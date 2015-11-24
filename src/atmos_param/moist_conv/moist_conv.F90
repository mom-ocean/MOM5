
module moist_conv_mod

!-----------------------------------------------------------------------

 use           mpp_mod, only : mpp_pe,             &
                               mpp_root_pe,        &
                               stdlog
use   time_manager_mod, only : time_type
 use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
use  sat_vapor_pres_mod, ONLY: lookup_es_des, compute_qs, descomp
use mpp_mod,             only: input_nml_file
use             fms_mod, ONLY:  error_mesg, file_exist, open_namelist_file,  &
                                check_nml_error, close_file,        &
                                FATAL, WARNING, NOTE, mpp_pe, mpp_root_pe, &
                                write_version_number, stdlog
use       constants_mod, ONLY: HLv, HLs, cp_air, grav, rdgas, rvgas

use           fms_mod, only : write_version_number, ERROR_MESG, FATAL
use field_manager_mod, only : MODEL_ATMOS
use tracer_manager_mod, only : get_tracer_index,   &
                               get_number_tracers, &
                               get_tracer_names,   &
                               get_tracer_indices, &
                               query_method,       &
                               NO_TRACER

implicit none
private

!------- interfaces in this module ------------

public :: moist_conv, moist_conv_Init, moist_conv_end

!-----------------------------------------------------------------------
!---- namelist ----

 real :: HC   = 1.00
 real :: beta = 0.0
 real :: TOLmin=.02, TOLmax=.10
 integer :: ITSMOD=30
 logical :: do_simple =.false.

!----- note beta is the fraction of convective condensation that is
!----- detrained into a stratiform cloud

 namelist /moist_conv_nml/  HC, beta, TOLmin, TOLmax, ITSMOD, do_simple

!-----------------------------------------------------------------------
!---- VERSION NUMBER -----

 character(len=128) :: version = '$Id: moist_conv.F90,v 19.0 2012/01/06 20:10:09 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
 logical            :: module_is_initialized = .false.

!---------- initialize constants used by this module -------------------

 real, parameter :: d622 = rdgas/rvgas
 real, parameter :: d378 = 1.0-d622
 real, parameter :: grav_inv = 1.0/grav
 real, parameter :: rocp = rdgas/cp_air

real :: missing_value = -999.
integer :: nsphum, nql, nqi, nqa   ! tracer indices for stratiform clouds

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_qldt_conv,   id_qidt_conv,   id_qadt_conv, &
           id_ql_conv_col, id_qi_conv_col, id_qa_conv_col,&
           id_q_conv_col, id_t_conv_col

character(len=3) :: mod_name = 'mca'

logical :: do_mca_tracer = .false.
integer :: num_mca_tracers = 0
integer               :: num_tracers
integer, allocatable, dimension(:) :: id_tracer_conv, id_tracer_conv_col

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS

!#######################################################################

 subroutine moist_conv ( Tin, Qin, Pfull, Phalf, coldT,        & ! required
                         Tdel, Qdel, Rain, Snow,               & ! required
                         dtinv, Time, is, js, tracers, qtrmca, & ! required
                         Lbot, mask, Conv,                     & ! optional
                         ql, qi, cf, qldel, qidel, cfdel)        ! optional

!-----------------------------------------------------------------------
!
!                       MOIST CONVECTIVE ADJUSTMENT
!
!-----------------------------------------------------------------------
!
!   INPUT:   Tin     temperature at full model levels
!            Qin     specific humidity of water vapor at full
!                      model levels
!            Pfull   pressure at full model levels
!            Phalf   pressure at half model levels
!            coldT   Should MCA produce snow in this column?
!
!   OUTPUT:  Tdel    temperature adjustment at full model levels (deg k)
!            Qdel    specific humidity adjustment of water vapor at
!                       full model levels
!            Rain    liquid precipitiation (in Kg m-2)
!            Snow    ice phase precipitation (kg m-2)
!  OPTIONAL
!
!   INPUT:   Lbot    integer index of the lowest model level,
!                      Lbot is always <= size(Tin,3)
!              ql    liquid water condensate
!              qi    ice condensate
!              cf    stratiform cloud fraction (used only when
!                    operating with stratiform cloud scheme) (fraction)
!
!  OUTPUT:   Conv    logical flag; TRUE then moist convective
!                       adjustment was performed at that model level.
!            cfdel   change in stratiform cloud fraction (fraction)
!            qldel   change in liquid water condensate due to
!                    convective detrainment (kg condensate /kg air)
!            qidel   change in ice condensate due to
!                    convective detrainment (kg condensate /kg air)
!
!-----------------------------------------------------------------------
!----------------------PUBLIC INTERFACE ARRAYS--------------------------
    real, intent(INOUT), dimension(:,:,:)           :: Tin, Qin
    real, intent(IN) ,   dimension(:,:,:)           :: Pfull, Phalf 
 logical, intent(IN) ,   dimension(:,:)             :: coldT
    real, intent(OUT),   dimension(:,:,:)           :: Tdel, Qdel
    real, intent(OUT),   dimension(:,:)             :: Rain, Snow
    real, intent(IN)                                :: dtinv
type(time_type), intent(in)                         :: Time
integer, intent(IN)                                :: is, js
    real, dimension(:,:,:,:), intent(in)            :: tracers
    real, dimension(:,:,:,:), intent(out)           :: qtrmca
 integer, intent(IN) ,   dimension(:,:),   optional :: Lbot
    real, intent(IN) ,   dimension(:,:,:), optional :: mask
 logical, intent(OUT),   dimension(:,:,:), optional :: Conv
    real, intent(INOUT), dimension(:,:,:), optional :: ql, qi, cf
    real, intent(OUT),   dimension(:,:,:), optional :: qldel, qidel, cfdel
         
!-----------------------------------------------------------------------
!----------------------PRIVATE (LOCAL) ARRAYS---------------------------
! logical, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) :: DO_ADJUST
!    real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) ::  &
!-----------------------------------------------------------------------
!----------------------PRIVATE (LOCAL) ARRAYS---------------------------
integer, dimension(size(Tin,1),size(Tin,2)) :: ISMVF
integer, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) :: IVF

real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) ::   &
  Qdif,Temp,Qmix,Esat,Qsat,Test1,Test2, &
   Esdiff_v

real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)-1) ::   &
  Thalf,DelPoP,Esm,Esd,ALRM

real, dimension(size(Tin,3)) :: C,Ta,Qa
real, dimension(size(Tin,1),size(Tin,2)) :: HL

integer :: i,j,k,kk,KX,ITER,MXLEV,MXLEV1,kstart,KTOP,KBOT,KBOTM1
real    :: ALTOL,Sum0,Sum1,Sum2,EsDiff,EsVal,Thaf,Pdelta

logical  :: cloud_tracers_present
real, dimension(size(Phalf,1),size(Phalf,2),size(Phalf,3)) :: pmass
real, dimension(size(Phalf,1),size(Phalf,2)) :: tempdiag
integer  :: tr, num_cld_tracers
logical :: used
!-----------------------------------------------------------------------

      if (.not. module_is_initialized) call ERROR_MESG( 'MCA',  &
                                 'moist_conv_init has not been called', FATAL )

      num_cld_tracers = count( (/present(ql),present(qi),present(cf),present(qldel),present(qidel),present(cfdel)/) )
      if(num_cld_tracers == 0) then
        cloud_tracers_present = .false.
      else if(num_cld_tracers == 6) then
        cloud_tracers_present = .true.
      else
        call error_mesg('moist_conv','Either all or none of the cloud tracers and their tendencies must be present',FATAL)
      endif

        do k=1,size(Tin,3)
          pmass(:,:,k) = (Phalf(:,:,k+1)-Phalf(:,:,k))/GRAV
        end do

      KX=size(Tin,3)

!------ compute Proper HL
      if(do_simple) then
            HL = HLv
      else
        WHERE (coldT)
              HL = HLs
        ELSEWHERE
              HL = HLv
        END WHERE
      endif

!------ convert spec hum to mixing ratio ------
      Temp(:,:,:)=Tin(:,:,:)
      Qmix(:,:,:)=Qin(:,:,:)

      do k=1,KX-1
         DelPoP(:,:,k)=(Pfull(:,:,k+1)-Pfull(:,:,k))/Phalf(:,:,k+1)
      enddo

!-------------SATURATION VAPOR PRESSURE FROM ETABL----------------------
!  compute qs; also return dqsdT

      call compute_qs (Temp, Pfull, Qsat, dqsdT=Esdiff_v, hc = hc, &
                                                            esat=Esat)
      Qdif(:,:,:)=Max(0.0,Qmix(:,:,:)-Qsat(:,:,:))

!-----------------------------------------------------------------------
!                  MOIST CONVECTIVE ADJUSTMENT
!-----------------------------------------------------------------------

!  *** Set initial tolerance ***

           ALTOL = TOLmin

      do k=1,KX-1
         Thalf(:,:,k)=0.50*(Temp(:,:,k)+Temp(:,:,k+1))
      enddo

      call lookup_es_des (Thalf, Esm, Esd)

      do k=1,KX-1
         ALRM(:,:,k)=rocp*DelPoP(:,:,k)*Thalf(:,:,k)  &
         *(Phalf(:,:,k+1)+d622*HL(:,:)*Esm(:,:,k)/Thalf(:,:,k)/rdgas)  &
         /(Phalf(:,:,k+1)+d622*HL(:,:)*Esd(:,:,k)/cp_air)
      enddo

      IVF  (:,:,KX)=0
      Test1(:,:,KX)=0.0
      Test2(:,:,KX)=0.0

      do k=1,KX-1
         Test1(:,:,k)=Temp(:,:,k+1)-Temp(:,:,k)
         Test2(:,:,k)=ALRM(:,:,k)+ALTOL-Test1(:,:,k)
      enddo

!!!!! Test1(:,:,:)=0.0-Qdif(:,:,:)
      Test1(:,:,:)=(0.0-Qdif(:,:,:))*Qsat(:,:,:)

!-------IVF=1 in unstable layers where both levels are saturated--------

      do k=1,KX-1
         where (Test1(:,:,k) < 0.0 .and. Test1(:,:,k+1) < 0.0 .and. &
                Test2(:,:,k) < 0.0)
                         IVF(:,:,k)=1
         elsewhere
                         IVF(:,:,k)=0
         endwhere
      enddo

!  ------ Set convection flag (for optional output only) --------

      if (Present(Conv)) then
         Conv(:,:,1)=(IVF(:,:,1) == 1)
         do k=1,KX-1
            Conv(:,:,k+1)=(IVF(:,:,k) == 1 .or. IVF(:,:,k+1) == 1)
         enddo
      endif

!  ----- Set counter for each column -----

         ISMVF(:,:)=0
      do k=1,KX-1
         ISMVF(:,:)=ISMVF(:,:)+IVF(:,:,k)
      enddo

!-----------------------------------------------------------------------
!---------------LOOP OVER EACH VERTICAL COLUMN--------------------------
                       do j=1,size(Tin,2)
           OUTER_LOOP: do i=1,size(Tin,1)
!-----------------------------------------------------------------------
      if (ISMVF(i,j) == 0)  CYCLE

      if (Present(Lbot)) then
         MXLEV=Lbot(i,j)
      else
         MXLEV=KX
      endif
      MXLEV1=MXLEV-1

!  *** Re-set initial tolerance ***
           ALTOL = TOLmin

!----------(return here after increasing tolerance)--------------------
1450  CONTINUE

!--------------Iterations at the same tolerance-------------------------
                     do 1740 ITER=1,ITSMOD
!-----------------------------------------------------------------------
      kstart=1
 1500 CONTINUE
!-------------TEST TO DETERMINE UNSTABLE LAYER BLOCKS-------------------
!-------Find top (KTOP) and bottom (KBOT) of unstable layers------------
      do k=kstart,MXLEV1
          if (IVF(i,j,k) == 1)  GO TO 1505
      enddo
      CYCLE OUTER_LOOP
1505  KTOP=k

      do k=KTOP,MXLEV1
          if (IVF(i,j,k+1) == 0) then
             KBOT=k+1
             GO TO 1510
          endif
      enddo
      KBOT=MXLEV
1510  CONTINUE
!-----------------------------------------------------------------------

      KBOTM1=KBOT - 1
      Sum1=0.0
      Sum2=0.0
!-----------------------------------------------------------------------
                      do 1630 k=KTOP,KBOT
!-----------------------------------------------------------------------
      if(do_simple) then
        call DEsComp (Temp(i,j,k),EsDiff)
        C(k)=d622*HC*EsDiff/Pfull(i,j,k)
      else
        C(k)=Pfull(i,j,k)-d378*Esat(i,j,k)
        if (C(k) <= 0.0) then
          C(k)=0.0
        else
          C(k) = esdiff_v(i,j,k)
        endif
      endif  

      Sum0=0.0
      if (k == KBOT) GO TO 1625
      kk=k
1620  if (kk > KBOTM1) GO TO 1625
      Sum0=Sum0+ALRM(i,j,kk)
      kk=kk+1
      GO TO 1620

1625  CONTINUE
      Pdelta=Phalf(i,j,k+1)-Phalf(i,j,k)
      Sum1=Sum1 + Pdelta*((cp_air+HL(i,j)*C(k))*(Temp(i,j,k)+Sum0)+  &
                           HL(i,j)*(Qmix(i,j,k)-Qsat(i,j,k)))
      Sum2=Sum2 + Pdelta*(cp_air+HL(i,j)*C(k))
!-----------------------------------------------------------------------
1630                   CONTINUE
!-----------------------------------------------------------------------
      Ta(KBOT)=Sum1/Sum2
      k=KTOP
1645  if (k > KBOTM1) GO TO 1641
      Sum0=0.0
      kk=k
1640  if (kk > KBOTM1) GO TO 1642
      Sum0=Sum0+ALRM(i,j,kk)
      kk=kk+1
      GO TO 1640
1642  Ta(k)=Ta(KBOT)-Sum0
      k=k+1
      GO TO 1645

!---------UPDATE T,R,ES,Esm,Esd & Qsat FOR THE ADJUSTED POINTS----------

1641  do k=KTOP,KBOT
        Qa(k)=Qsat(i,j,k)+C(k)*(Ta(k)-Temp(i,j,k))
        Temp(i,j,k)=Ta(k)
        Qmix(i,j,k)=Qa(k)
        call compute_qs ( Temp(i,j,k), Pfull(i,j,k), Qsat(i,j,k), &
                     dqsdT = Esdiff_v(i,j,k), hc =hc, esat= Esat(i,j,k))
        Qdif(i,j,k)=Max(0.0,Qmix(i,j,k)-Qsat(i,j,k))
      enddo

      do k=KTOP,KBOTM1
        Thaf=0.50*(Temp(i,j,k)+Temp(i,j,k+1))
!DIR$ INLINE
        call lookup_es_des (Thaf, EsVal, Esdiff)
!DIR$ NOINLINE
        Esm (i,j,k)=HC*EsVal
        Esd (i,j,k)=HC*EsDiff
        ALRM(i,j,k)=rocp*DelPoP(i,j,k)*Thaf*  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esm(i,j,k)/Thaf/rdgas)/  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esd(i,j,k)/cp_air)
      enddo

!------------Is this the bottom of the current column ???---------------
      kstart=KBOT+1
      if (kstart <= MXLEV1) GO TO 1500
!-----------------------------------------------------------------------
                if (ITER == ITSMOD) GO TO 1740
!-----------------------------------------------------------------------

      do k=1,MXLEV1
        Thaf=0.50*(Temp(i,j,k)+Temp(i,j,k+1))
!DIR$ INLINE
        call lookup_es_des (Thaf, EsVal, EsDiff)
!DIR$ NOINLINE
        Esm (i,j,k)=HC*EsVal
        Esd (i,j,k)=HC*EsDiff
        ALRM(i,j,k)=rocp*DelPoP(i,j,k)*Thaf*  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esm(i,j,k)/Thaf/rdgas)/  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esd(i,j,k)/cp_air)
      enddo

      do k=1,MXLEV1
        IVF(i,j,k)=0
!!!!    if (Qdif(i,j,k) > 0.0 .and. Qdif(i,j,k+1) > 0.0 .and.  &
        if (Qdif(i,j,k)*Qsat(i,j,k) > 0.0 .and.     &
             Qdif(i,j,k+1)*Qsat(i,j,k+1) > 0.0 .and.  &
              (Temp(i,j,k+1)-Temp(i,j,k)) > (ALRM(i,j,k)+ALTOL)) then
                       IVF(i,j,k) = 1
        endif
      enddo

!   ------ reset optional convection flag ------

      if (Present(Conv)) then
         Conv(i,j,1)=(IVF(i,j,1) == 1)
         do k=1,MXLEV1
            Conv(i,j,k+1)=(IVF(i,j,k) == 1 .or. IVF(i,j,k+1) == 1)
         enddo
      endif

!   ------ Are all layers sufficiently stable ??? ------

              ISMVF(i,j)=0
           do k=1,MXLEV1
              ISMVF(i,j)=ISMVF(i,j)+IVF(i,j,k)
           enddo
              if (ISMVF(i,j) == 0) CYCLE OUTER_LOOP

!-----------------------------------------------------------------------
1740                       CONTINUE
!-----------------------------------------------------------------------

!---------Maximum iterations reached: Increase tolerance (ALTOL)--------
      ALTOL=2.0*ALTOL
!del  WRITE (*,9902) I,ALTOL
      call error_mesg ('moist_conv', 'Tolerence (ALTOL) doubled', NOTE)
      if (ALTOL <= TOLmax)  GO TO 1450

!     WRITE (*,9903)
!     WRITE (*,9904) (k,Temp(i,j,k),Qmix(i,j,k),Qsat(i,j,k),  &
!                       Qdif(i,j,k),ALRM(i,j,k),k=1,MXLEV1)
!     WRITE (*,9904) (k,Temp(i,j,k),Qmix(i,j,k),Qsat(i,j,k),  &
!                       Qdif(i,j,k)            ,k=MXLEV,MXLEV)

   call error_mesg ('moist_conv', 'maximum iterations reached', WARNING)
!-----------------------------------------------------------------------
                            enddo OUTER_LOOP
                            enddo
!-----------------------------------------------------------------------
!---------------------- END OF i,j LOOP --------------------------------
!-----------------------------------------------------------------------

!---- call Convective Detrainment subroutine -----

     if (cloud_tracers_present) then

          !reset quantities
          cfdel = 0.
          qldel = 0.
          qidel = 0.
     
          CALL CONV_DETR(Qmix,Qin,Phalf,Temp,cf,coldT,cfdel,qldel,qidel)
          
     endif

!----- compute adjustments to temp and spec hum ----

      Tdel(:,:,:)=Temp(:,:,:)-Tin(:,:,:)
      Qdel(:,:,:)=Qmix(:,:,:)-Qin(:,:,:)

!----- integrate precip -----

      Rain(:,:)=0.0
      Snow(:,:)=0.0
   do k =1,KX

     if(do_simple) then
       Rain(:,:)=Rain(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               Qdel(:,:,k)*grav_inv
       if (cloud_tracers_present) then
          Rain(:,:)=Rain(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                                qldel(:,:,k)*grav_inv
       endif
     else
       WHERE(coldT(:,:)) 
         Snow(:,:)=Snow(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               Qdel(:,:,k)*grav_inv
       ELSEWHERE
         Rain(:,:)=Rain(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               Qdel(:,:,k)*grav_inv
       END WHERE

      !subtract off detrained condensate from surface precip
       if (cloud_tracers_present) then
         WHERE(coldT(:,:)) 
           Snow(:,:)=Snow(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                                qidel(:,:,k)*grav_inv
         ELSEWHERE
           Rain(:,:)=Rain(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                                qldel(:,:,k)*grav_inv
         END WHERE      
       end if
     endif

   enddo
      Rain(:,:)=Max(Rain(:,:),0.0)
      Snow(:,:)=Max(Snow(:,:),0.0)
!-----------------------------------------------------------------------
!-----------------   PRINT FORMATS   -----------------------------------

 9902 FORMAT(    ' *** ALTOL DOUBLED IN CONVAD AT I=',  &
                 I5,' ,ALTOL=', F10.4 )
 9903 FORMAT(/,' *** DIVERGENCE IN MOIST CONVECTIVE ADJUSTMENT ',/,  &
           4X,'K',14X,'T',14X,'R',13X,'Qsat',14X,'Qdif',12X,'ALRM',/)
 9904 FORMAT (I5,5E15.7)
!-----------------------------------------------------------------------


!------- update input values and compute tendency -------
                    
      Tin=Tin+Tdel;    Qin=Qin+Qdel
      
      Tdel=Tdel*dtinv; Qdel=Qdel*dtinv
      Rain=Rain*dtinv; Snow=Snow*dtinv
!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (cloud_tracers_present) then
         ql(:,:,:)=ql(:,:,:)+qldel(:,:,:)
         qi(:,:,:)=qi(:,:,:)+qidel(:,:,:)
         cf(:,:,:)=cf(:,:,:)+cfdel(:,:,:)
 
         qldel(:,:,:)=qldel(:,:,:)*dtinv
         qidel(:,:,:)=qidel(:,:,:)*dtinv
         cfdel(:,:,:)=cfdel(:,:,:)*dtinv
      endif   
      
!---------------------------------------------------------------------
!   define the effect of moist convective adjustment on the tracer 
!   fields. code to do so does not currently exist.
!---------------------------------------------------------------------
      qtrmca = 0.
 

!------- diagnostics for dt/dt_ras -------
      if ( id_tdt_conv > 0 ) then
        used = send_data ( id_tdt_conv, Tdel, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_ras -------
      if ( id_qdt_conv > 0 ) then
        used = send_data ( id_qdt_conv, Qdel, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_ras -------
      if ( id_prec_conv > 0 ) then
        used = send_data ( id_prec_conv, Rain+Snow, Time, is, js )
      endif
!------- diagnostics for snow_ras -------
      if ( id_snow_conv > 0 ) then
        used = send_data ( id_snow_conv, Snow, Time, is, js )
      endif

!------- diagnostics for water vapor path tendency ----------
      if ( id_q_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + Qdel(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_q_conv_col, tempdiag, Time, is, js )
      end if
   
!------- diagnostics for dry static energy tendency ---------
      if ( id_t_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + Tdel(:,:,k)*cp_air*pmass(:,:,k)
        end do
        used = send_data ( id_t_conv_col, tempdiag, Time, is, js )
      end if
   
   !------- stratiform cloud tendencies from cumulus convection ------------
   if (cloud_tracers_present) then

      !------- diagnostics for dql/dt from RAS or donner -------
      if ( id_qldt_conv > 0 ) then
        used = send_data ( id_qldt_conv, qldel(:,:,:), Time, is, js, 1, &
                           rmask=mask )
      endif
      
      !------- diagnostics for dqi/dt from RAS or donner -------
      if ( id_qidt_conv > 0 ) then
        used = send_data ( id_qidt_conv, qidel(:,:,:), Time, is, js, 1, &
                           rmask=mask )
      endif
      
      !------- diagnostics for dqa/dt from RAS or donner -------
      if ( id_qadt_conv > 0 ) then
        used = send_data ( id_qadt_conv, cfdel(:,:,:), Time, is, js, 1, &
                           rmask=mask )
      endif

      !------- diagnostics for liquid water path tendency ------
      if ( id_ql_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qldel(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_ql_conv_col, tempdiag, Time, is, js )
      end if
      
      !------- diagnostics for ice water path tendency ---------
      if ( id_qi_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qidel(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_qi_conv_col, tempdiag, Time, is, js )
      end if
      
      !---- diagnostics for column integrated cloud mass tendency ---
      if ( id_qa_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + cfdel(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_qa_conv_col, tempdiag, Time, is, js )
      end if
         
   end if ! if ( cloud_tracers_present )

   do tr = 1, num_mca_tracers
!------- diagnostics for dtracer/dt from RAS -------------
     if ( id_tracer_conv(tr) > 0 ) then
       used = send_data ( id_tracer_conv(tr), qtrmca(:,:,:,tr), Time, is, js, 1, &
                          rmask=mask )
     endif
 
!------- diagnostics for column tracer path tendency -----
     if ( id_tracer_conv_col(tr) > 0 ) then
       tempdiag(:,:)=0.
       do k=1,kx
         tempdiag(:,:) = tempdiag(:,:) + qtrmca(:,:,k,tr)*pmass(:,:,k)
       end do
       used = send_data ( id_tracer_conv_col(tr), tempdiag, Time, is, js )
     end if

 
   enddo

 end subroutine moist_conv

!#######################################################################
!#######################################################################

SUBROUTINE CONV_DETR(qvout,qvin,phalf,T,cf,coldT,cfdel,qldel,qidel)


IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a fraction of the water condensed
!      by the convection scheme and detrains in the top level
!      undergoing convective adjustment.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       VARIABLES
!
!   ------
!   INPUT:
!   ------
!
!       qvout    water vapor specific humidity AFTER adjustment
!                (kg vapor/kg air)
!       qvin     water  vapor specific humidity BEFORE adjustment
!                (kg vapor/kg air)
!       phalf    pressure at model half levels (Pascals)
!       T        Temperature (Kelvin)
!       cf       cloud fraction (fraction)
!       coldT    is condensation of ice nature?
!
!   -------------
!   INPUT/OUTPUT:
!   -------------
!
!       cfdel    Change in cloud fraction due to detrainment (fraction)
!       qldel    Increase in liquid water due to detrainment
!                (kg condensate/kg air)
!       qidel    Increase in ice due to detrainment
!                (kg condensate/kg air)
!
!   -------------------
!       INTERNAL VARIABLES:
!   -------------------
!
!       precipsource  accumulated source of precipitation
!                     (kg condensate /meter/ (seconds*squared))
!       ktop          integer of top level undergoing convection
!       accum         logical variable indicating whether or not to
!                     add precip
!       fT            fraction of condensate that is liquid
!       i,j,k         looping variables
!       IDIM,JDIM,KDIM dimensions of input arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:,:)  :: qvout,qvin,T,cf
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
LOGICAL,  INTENT (IN), DIMENSION(:,:)    :: coldT
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cfdel,qldel,qidel

!  Internal variables
!  ------------------

INTEGER                                  :: i,j,k,IDIM,JDIM,KDIM,ktop
LOGICAL                                  :: accum
REAL                                     :: precipsource

!
! Code
! ----

        ! reinitialize variables
        cfdel(:,:,:)   = 0.
        qidel(:,:,:)   = 0.
        qldel(:,:,:)   = 0.
        IDIM           = SIZE(qvout,1)
        JDIM           = SIZE(qvout,2)
        KDIM           = SIZE(qvout,3)

        !---loop over grid columns----!
        DO i = 1, IDIM
        DO j = 1, JDIM

             !reset variables
             precipsource     = 0.
             accum            = .FALSE.

             DO k = 1, KDIM

                 !begin new convective event
                 IF ((qvout(i,j,k) .ne. qvin(i,j,k)) .and. &
                     (.NOT. accum)) THEN
                     ktop  = k
                     accum = .TRUE.
                 END IF

                 !if convective event is over compute detrainment
                 IF ( (accum) .and. (qvout(i,j,k) .eq. qvin(i,j,k)) &
                      .and. (precipsource .gt. 0.)) THEN                      
                      if (coldT(i,j)) then
                      qidel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                      else
                      qldel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                      end if
                      cfdel(i,j,ktop) = MAX(0.,HC-cf(i,j,ktop))
                      accum        = .FALSE.
                      precipsource = 0.
                 END IF

                 !accumulate precip
                 IF (accum) THEN
                      precipsource = precipsource + &
                                   ( qvin(i,j,k)  -qvout(i,j,k))* &
                                   (phalf(i,j,k+1)-phalf(i,j,k))
                 END IF

             END DO    !---end k loop over vertical column

            !---clear any remaining precip
            IF ( (precipsource .gt. 0.) .and. (accum) ) THEN
                 if (coldT(i,j)) then
                 qidel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                 else
                 qldel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                 end if
            END IF

        END DO    !---end j loop
        END DO    !---end i loop

END SUBROUTINE CONV_DETR

!#######################################################################

subroutine moist_conv_init (axes, Time, tracers_in_mca)

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 logical, dimension(:), intent(in), optional :: tracers_in_mca

!-----------------------------------------------------------------------
      
 integer :: unit, io, ierr, logunit
 integer :: nn, tr
 character(len=128) :: diagname, diaglname, tendunits, name, units


!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=moist_conv_nml, iostat=io)
    ierr = check_nml_error(io,"moist_conv_nml")
#else
    if (file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
            read  (unit, nml=moist_conv_nml, iostat=io, end=10)
            ierr = check_nml_error (io,'moist_conv_nml')
        enddo
 10     call close_file (unit)
    endif
#endif

!---------- output namelist --------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit,nml=moist_conv_nml)
      endif

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from moist conv adj',     'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from moist conv adj',   'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from moist conv adj',       'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from moist conv adj',       'kg/m2/s' )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from moist conv adj','kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from moist conv adj','W/m2' )


!---------------------------------------------------------------------
! --- Find the tracer indices 
!---------------------------------------------------------------------
  call get_number_tracers (MODEL_ATMOS, num_prog= num_tracers)
  if ( num_tracers .gt. 0 ) then
  else
    call error_mesg('moist_conv_init', 'No atmospheric tracers found', FATAL)
  endif
    ! get tracer indices for stratiform cloud variables
      nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )

     
!----------------------------------------------------------------------
!    determine how many tracers are to be transported by moist_conv_mod.
!----------------------------------------------------------------------
      num_mca_tracers = count(tracers_in_mca)
      if (num_mca_tracers > 0) then
        do_mca_tracer = .true.
      else
        do_mca_tracer = .false.
      endif

!---------------------------------------------------------------------
!    allocate the arrays to hold the diagnostics for the moist_conv 
!    tracers.
!---------------------------------------------------------------------
      allocate(id_tracer_conv    (num_mca_tracers)) ; id_tracer_conv = 0
      allocate(id_tracer_conv_col(num_mca_tracers)) ; id_tracer_conv_col = 0
      nn = 1
      do tr = 1,num_tracers
        if (tracers_in_mca(tr)) then
          call get_tracer_names(MODEL_ATMOS, tr, name=name, units=units)
 
!----------------------------------------------------------------------
!    for the column tendencies, the name for the diagnostic will be 
!    the name of the tracer followed by 'dt_MCA'. the longname will be 
!    the name of the tracer followed by ' tendency from MCA'. units are
!    the supplied units of the tracer divided by seconds.
!----------------------------------------------------------------------
      diagname = trim(name)//'dt_MCA'
      diaglname = trim(name)//' tendency from MCA'
      tendunits = trim(units)//'/s'
      id_tracer_conv(nn) = register_diag_field ( mod_name, &
                             trim(diagname), axes(1:3), Time, &
                             trim(diaglname), trim(tendunits),  &
                             missing_value=missing_value        )

!----------------------------------------------------------------------
!    for the column integral  tendencies, the name for the diagnostic 
!    will be the name of the tracer followed by 'dt_MCA_col'. the long-
!    name will be the name of the tracer followed by ' path tendency 
!    from MCA'. units are the supplied units of the tracer multiplied
!    by m**2 /kg divided by seconds.
!----------------------------------------------------------------------
      diagname = trim(name)//'dt_MCA_col'
      diaglname = trim(name)//' path tendency from MCA'
      tendunits = trim(units)//'m2/kg/s'
      id_tracer_conv_col(nn) = register_diag_field ( mod_name, &
                                 trim(diagname), axes(1:2), Time, &
                                 trim(diaglname), trim(tendunits),  &
                                 missing_value=missing_value)
      nn = nn + 1
     endif
    end do




      module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine moist_conv_init


!#######################################################################
subroutine moist_conv_end

integer :: log_unit

if(.not.module_is_initialized) then
  return
else
  module_is_initialized = .FALSE.
endif

log_unit = stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (log_unit,'(/,(a))') 'Exiting moist_conv.'
endif

end subroutine moist_conv_end

!#######################################################################

end module moist_conv_mod


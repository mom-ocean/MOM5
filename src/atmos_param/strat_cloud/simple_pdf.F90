MODULE simple_pdf_mod

!MNS NOTE: THIS IS EXPERIMENTAL...

use fms_mod,                   only: write_version_number
use beta_dist_mod,             only: incomplete_beta, beta_dist_init,  &
                                     beta_dist_end
use strat_cloud_utilities_mod, only: strat_cloud_utilities_init, &
                                     diag_id_type, diag_pt_type

implicit none
private

!-----------------------------------------------------------------------
!----interfaces--------------------------------------------------------
public  simple_pdf, simple_pdf_init, simple_pdf_end

!----------------------------------------------------------------------
!----version number----------------------------------------------------
Character(len=128) :: Version = '$Id: simple_pdf.F90,v 20.0 2013/12/13 23:22:07 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'



logical            :: module_is_initialized = .false.



CONTAINS



!#########################################################################

SUBROUTINE simple_pdf_init

!-----------------------------------------------------------------------
      if (module_is_initialized) return

!-----------------------------------------------------------------------
!    write version number to output file.
!-----------------------------------------------------------------------
      call write_version_number (version, tagname)

!-----------------------------------------------------------------------
!    make sure needed modules have been initialized.
!-----------------------------------------------------------------------
      call strat_cloud_utilities_init
      call beta_dist_init

!-----------------------------------------------------------------------
!    mark this module initialized.
!-----------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------

END SUBROUTINE simple_pdf_init


!#########################################################################

SUBROUTINE simple_pdf (j, idim, jdim, kdim, qmin, qa, qtot, qs, gamma, &
                       qthalfwidth, betaP, inv_dtcloud, SA_0, n_diag_4d, &
                       diag_4d, diag_id, diag_pt, SA, qa_upd)    

!-------------------------------------------------------------------------
INTEGER,                               INTENT(IN)    :: j, idim, jdim, kdim
REAL,                                  INTENT(IN)    :: qmin, inv_dtcloud
REAL, dimension(idim,kdim),            INTENT(IN)    :: qtot, qa, qs, &
                                                        gamma, SA_0
INTEGER,                               INTENT(IN)    :: betaP
REAL,                                  INTENT(IN)    :: qthalfwidth
INTEGER,                               INTENT(IN)    :: n_diag_4d
REAL, dimension(idim, jdim, kdim, 0:n_diag_4d),   &
                                       INTENT(INOUT) :: diag_4d
TYPE(diag_id_type),                    intent(in)    :: diag_id
TYPE(diag_pt_type),                    intent(inout) :: diag_pt
REAL, dimension(idim, kdim),           INTENT(OUT)   :: SA, qa_upd

!------------------------------------------------------------------------
!local variables
      REAL, dimension(idim, kdim)   :: qta3, qtqsa3, qcg, qag, qa1, qa0, &
                                       qabar
      real, dimension(idim)         :: qtbar, deltaQ, qtmin, qs_norm, &
                                       qagtmp, qcgtmp, qvgtmp
      real                          :: icbp, icbp1, pnorm
      INTEGER                       :: k, id

!------------------------------------------------------------------------
!    compute pdf cloud fraction and condensate. Note that the SYMMETRIC 
!    beta distribution is used here. Initialize grid-box mean values of 
!    cloud fraction (qag) and cloud condensate(qcg).
!------------------------------------------------------------------------
      qta3(:,:) = max(qmin, qtot)
      qtqsa3(:,:) = qta3(:,:) - qs
      qcg(:,:)  = 0.
      qag(:,:)  = 0.
      DO k=1,kdim
        
!------------------------------------------------------------------------
!    calculate normalized vertical level. 0. = top of gridbox, 1. = bottom 
!    of gridbox
!------------------------------------------------------------------------
        pnorm =  0.5 
        
!------------------------------------------------------------------------
!    First step is to calculating the minimum (qtmin) of the total water 
!    distribution and the width of the qt distribution (deltaQ). For 
!    diagnostic variance this is set to (1.-qthalfwidth)*qtbar and 
!    2*qthalfwidth*qtbar, respectively, where qtbar is the mean total water
!    in the grid box.        
!------------------------------------------------------------------------
        qtbar = qta3(:,k)
        qtbar  = max(qmin ,qtbar)
        deltaQ  = 2.*qthalfwidth*qtbar
        qtmin = (1. - qthalfwidth)*qtbar 
        
!------------------------------------------------------------------------
!    From this the variable normalized saturation specific humidity qs_norm
!    is  calculated.
!    qs_norm = (qs(Tl) - qtmin)/(qtmax-qtmin)
!
!          = 0.5  - (qtbar - qs(Tl))/deltaQ
!
!    Note that if qs_norm > 1., the grid box is fully clear. If qs_norm 
!    < 0., the grid box is fully cloudy.
!------------------------------------------------------------------------
        qs_norm = qtqsa3(:,k)
        qs_norm = 0.5 - (qs_norm/deltaQ)
             
!------------------------------------------------------------------------
!    Calculation of cloud fraction (qagtmp), cloud condensate (qcgtmp), and
!    water vapor in clear air part of the grid box (qvgtmp).
!
!    Formulas (from Tompkins, and personal derivations):
!
!    Define icbp  = incomplete_beta(qs_norm,p,q)
!           icbp1 = incomplete_beta(qs_norm,p+1,q)
!
!    qagtmp = 1. - icbp
!
!    qcgtmp = aThermo * {  (qtbar-qtmin)*(1.-icbp1) - 
!                       qs_norm*deltaQ*(1.-icbp ) }
!
!
!    qvgtmp = qtmin + (p/(p+q))*(icbp1/icbp)*deltaQ
!
!  
!    where aThermo = 1./(1.+(L/cp)*dqsdT)
!
!    note that in the qvg formula below the factor of 0.5 is equal 
!    to (p/(p+q)).
!
!------------------------------------------------------------------------
        do id = 1,idim
          if (qs_norm(id).le.1.) then
            icbp = incomplete_beta (max(0., qs_norm(id)), &
                                    p = betaP, q = betaP)
            icbp1 = incomplete_beta (max(0., qs_norm(id)), &
                                     p = betaP + 1, q = betaP)
            qagtmp(id) = 1. - icbp
            qcgtmp(id) = (qtbar(id) - qtmin(id))*(1. - icbp1) -  &
                           qs_norm(id)*deltaQ(id)*(1. - icbp)    
            qcgtmp(id) = qcgtmp(id)/(1. + gamma(id,k))
            qvgtmp(id) = qtmin(id) + 0.5*(icbp1/max(icbp, qmin))*deltaQ(id)
             
!------------------------------------------------------------------------
!    bound very very small cloud fractions which may cause negative cloud 
!    condensates due to roundoff errors or similar errors in the beta table 
!    lookup.
!------------------------------------------------------------------------
            if ((qagtmp(id) .lt. 0.).or.(qcgtmp(id) .le. 0.)) then
              qagtmp(id) = 0.
              qcgtmp(id) = 0.
              qvgtmp(id) = qtbar(id)
            end if
          else             
            qagtmp(id) = 0.
            qcgtmp(id) = 0.
            qvgtmp(id) = qtbar(id)             
          end if
        end do
              
!------------------------------------------------------------------------
!    sum vertically. note special averaging of clear-sky water vapor
!    this is weighting clear-sky relative humidity by the clear-sky fraction
!------------------------------------------------------------------------
        qag(:,k) = qag(:,k) + qagtmp
        qcg(:,k) = qcg(:,k) + qcgtmp
      enddo ! ( k loop) 

!------------------------------------------------------------------------
!    compute grid-box average cloud fraction, cloud condensate and 
!    water vapor.  do adjustment of cloud fraction. set total tendency term
!    (SA)  and update cloud fraction (qa_upd).    
!------------------------------------------------------------------------
      qa0 = qa
      qa1 = qag
      SA  = SA_0 + qa1 - qa0
      qa_upd = qa1

!--------------------------------------------------------------------------
!    fill desired diagnostics.
!--------------------------------------------------------------------------
      if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
        diag_4d(:,j,:,diag_pt%qadt_lsform ) =    &
                                           max(qa1 - qa0, 0.)*inv_dtcloud 
      end if
      if (diag_id%qadt_lsdiss + diag_id%qa_lsdiss_col > 0) then
        diag_4d(:,j,:,diag_pt%qadt_lsdiss ) =    &
                                             max(qa0 - qa1, 0.)*inv_dtcloud
      end if

!------------------------------------------------------------------------

END SUBROUTINE simple_pdf


!#########################################################################

SUBROUTINE simple_pdf_end

      call beta_dist_end

      module_is_initialized = .false.

END SUBROUTINE simple_pdf_end


!#########################################################################



END MODULE simple_pdf_mod

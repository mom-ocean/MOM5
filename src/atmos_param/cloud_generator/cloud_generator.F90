module cloud_generator_mod

!   shared modules:
  use sat_vapor_pres_mod, only: compute_qs
  use constants_mod,      only: hlv, hls, cp_air, tfreeze, &
                                rvgas, rdgas
  use mpp_mod,            only: input_nml_file
  use fms_mod,            only: open_namelist_file, mpp_pe,       &
                                mpp_root_pe, stdlog,              &
                                write_version_number, file_exist, &
                                check_nml_error, error_mesg,      &
                                FATAL, close_file
  use random_numbers_mod, only: randomNumberStream,        &
                                getRandomNumbers
  use beta_dist_mod,      only: beta_dist_init, beta_dist_end, &
                                incomplete_beta, beta_deviate
!--------------------------------------------------------------------
  !
  ! Given a profile of cloud fraction, produce a set of columns indicating 
  !   the presence or absence of cloud consistent with one of four overlap  
  !   assumptions: random, maximum, maximum-random, and one allowing 
  !   the rank correlation of the variable to be specified between layers.
  ! The module uses a random number generation module which can be used to wrap
  !   an arbitrary random number generator. The module defines a type that keeps 
  !   the state (the seed, a generator indicator, or whatever).  
  ! Each function takes cloud fraction as a function of height. You can either supply 
  !   a vector or a 3D array of dimensions (nX, ny, nLevels); in either case the 
  !   first element in the level dimension is the highest model layer. 
  ! Each function returns an array nSamples by nLevels long, where nLevels
  !   is determined by the size of the cloud fraction array. 
  ! The neighbor-to-neighbor correlation routine takes an additional 
  !   parameters described below. 
  !
!--------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloud_generator.F90,v 19.0 2012/01/06 20:02:09 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

  interface genRandomOverlapSamples
    module procedure genRandomOverlapSamples_1D, genRandomOverlapSamples_3D
  end interface ! genRandomOverlapSamples
  
  interface genMaximumOverlapSamples
    module procedure genMaximumOverlapSamples_1D, genMaximumOverlapSamples_3D
  end interface ! genMaximumOverlapSamples
  
  interface genMaxRanOverlapSamples
    module procedure genMaxRanOverlapSamples_1D, genMaxRanOverlapSamples_3D
  end interface ! genMaxRanOverlapSamples
  
  interface genWeightedOverlapSamples
    module procedure genWeightedOverlapSamples_1D, genWeightedOverlapSamples_3D
  end interface ! genWeightedOverlapSamples
  
  public :: cloud_generator_init, &
            cloud_generator_end,  &
            generate_stochastic_clouds, &
            do_cloud_generator,   &
            compute_overlap_weighting
  
!---------------------------------------------------------------------
!-------- namelist  ---------

  ! Minimum values for cloud fraction, water, ice contents
  !   Taken from cloud_rad. Perhaps these should be namelist parameters? 
  real, parameter :: qcmin = 1.E-10, qamin = 1.E-2
  
  ! Pressure scale height - for converting pressure differentials to height
  real, parameter :: pressureScaleHeight = 7.3 ! km
  
  ! Overlap parameter: 1 - Maximum, 2 - Random, 3 - Maximum/Random,
  !                    4 - Exponential
  integer         :: defaultOverlap = 2 
  real            :: overlapLengthScale = 2.0  ! km
  ! These control the option to pull cloud condensate from a symmetric  
  !    beta distribution with p = q = betaP, adjusting 
  !    the upper and lower bounds of the ditribution to match 
  !    cloud fraction and cloud condensate. 
  logical         :: do_inhomogeneous_clouds = .false. 
  integer         :: betaP = 5

  !The following apply for pdf cloud scheme
  !
  ! qthalfwidth defines the fraction of qtbar (mean total water in the
  ! grid box) that the maximum and minimum of the distribution differ 
  ! from qtbar. That is, total water at the sub-grid scale may take on 
  ! values anywhere between (1.-qthalfwidth)*qtbar and (1.+qthalfwidth)*
  ! qtbar
  !
  logical         :: do_pdf_clouds = .false.
  real            :: qthalfwidth = 0.1  
  
  ! The following apply for ppm vertical interpolation if done.
  !
  !      nsublevels      This is the number of sub-levels to be used
  !                      for sub-grid scale vertical structure to
  !                      clouds. If equal to 1, then no vertical
  !                      sub-grid scale structure is calculated.
  !
  !      kmap, kord      Quantities related to the PPM vertical inter-
  !                      polation calculation.
  !
  integer           :: nsublevels     = 1
  integer           :: kmap           = 1
  integer           :: kord           = 7
  
namelist /cloud_generator_nml/  defaultOverlap, overlapLengthScale, &
                                  do_inhomogeneous_clouds, betaP,     &
                                  do_pdf_clouds, qthalfwidth, nsublevels, &
                                  kmap,kord

!----------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.  ! module is initialized ?
logical :: cloud_generator_on    = .false.  ! is module being operated?

!----------------------------------------------------------------------

                              contains 
                              
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              
!######################################################################
subroutine cloud_generator_init


!---------------------------------------------------------------------
!    cloud_generator_init is the constructor for 
!    cloud_generator_mod.

!----------------------------------------------------------------------
!   local variables:
      integer   ::   unit, ierr, io, logunit

!--------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!
!--------------------------------------------------------------------

      if (.not. module_is_initialized) then
!---------------------------------------------------------------------
!    read namelist.         
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=cloud_generator_nml, iostat=io)
        ierr = check_nml_error(io,"cloud_generator_nml")
#else
        if (file_exist('input.nml')) then
          unit =  open_namelist_file ( )
          ierr=1; do while (ierr /= 0)
          read (unit, nml=cloud_generator_nml, iostat=io, end=10) 
          ierr = check_nml_error (io, 'cloud_generator_nml')
          enddo                       
10        call close_file (unit)      
        endif                         
#endif
!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
        call write_version_number(version, tagname)
        logunit = stdlog()
        if (mpp_pe() == mpp_root_pe() ) &
                   write (logunit, nml=cloud_generator_nml)
                   
!-----------------------------------------------------------------------
!    do_inhomogeneous_clouds and do_pdf_clouds cannot be 
!    simultaneously true
!-----------------------------------------------------------------------
        if (do_inhomogeneous_clouds .and. do_pdf_clouds) then
          call error_mesg ( 'cloud_generator_mod', &
           'do_inhomogeneous_clouds and do_pdf_clouds cannot'//&
           'be simultaneously true', FATAL)
        endif
!-----------------------------------------------------------------------
!    qthalfwidth must be greater than 0.
!-----------------------------------------------------------------------
        if (qthalfwidth .lt. 1.e-03) then
          call error_mesg ( 'cloud_generator_mod', &
           'qthalfwidth must be greater than 0.001', FATAL)
        endif

!---------------------------------------------------------------------
!    Initialize the beta distribution module if we're going to need it. 
!---------------------------------------------------------------------
        if(do_inhomogeneous_clouds .or. do_pdf_clouds) &
             call beta_dist_init

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
        module_is_initialized = .true.
        cloud_generator_on    = .true.
     end if

end subroutine cloud_generator_init
!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine generate_stochastic_clouds(streams, ql, qi, qa, qn, qni,    &
                                      overlap, pFull, pHalf, & 
                                      temperature, qv,&
                                      cld_thickness, &
                                      ql_stoch, qi_stoch, qa_stoch, &
                                      qn_stoch, qni_stoch)
!--------------------------------------------------------------------
!   intent(in) variables:
!
  type(randomNumberStream), &
           dimension(:, :),     intent(inout) :: streams
  ! Dimension nx, ny, nz
  real,    dimension(:, :, :),    intent( in) :: ql, qi, qa, qn, qni
  integer,                     optional, &
                                  intent( in) :: overlap
  real,    dimension(:, :, :), optional, &
                                 intent( in)  :: pFull, temperature, qv
  ! Dimension nx, ny, nz+1
  real,    dimension(:, :, :), optional, &
                                 intent( in)  :: pHalf                  
  ! Dimension nx, ny, nz, nCol = nBands
  integer, dimension(:, :, :, :), intent(out) :: cld_thickness 
  real,    dimension(:, :, :, :), intent(out) :: ql_stoch, qi_stoch, &
                                                 qa_stoch, qn_stoch, &
                                                 qni_stoch
  ! ---------------------------------------------------------
  ! Local variables
  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3)) :: qa_local, ql_local, qi_local
  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3)) :: heightDifference, &
                                           overlapWeighting  ! 1 for max, 0 for random
                         
  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3), &
                     size(ql_stoch, 4)) :: pdfRank ! continuously from 0 to 1
                                  
  ! These arrays could be declared allocatable and used only when 
  !    do_inhomogeneous_clouds OR do_pdf_clouds is true. 
  real,    dimension(size(ql_stoch, 1), &  ! Quantities for estimating condensate variability
                     size(ql_stoch, 2), &  !   from a beta distribution
                     size(ql_stoch, 3)) :: aThermo, qlqcRatio, qs_norm, &
                                           deltaQ, qs

  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3), &
                     size(ql_stoch, 4)) :: qc_stoch !stochastic condensate
  real,    dimension(4,                 &
                     size(ql_stoch,1),  &
                     size(ql_stoch,2),  &
                     size(ql_stoch,3)) :: qta4,qtqsa4 !used for ppm
  real,    dimension(size(ql_stoch,1),  &
                     size(ql_stoch,2),  &
                     size(ql_stoch,3)) :: delp !used for ppm
  real,    dimension(size(ql_stoch,1),  &
                     size(ql_stoch,2),  &
                     size(ql_stoch,4)) :: qctmp, qatmp !temporary summing 
                                                       !variables
                                   
  real    :: pnorm   !fraction of level- used for ppm interpolation    
  real    :: tmpr  
  integer :: nLev, nCol
  integer :: overlapToUse
  integer :: i, j, k, n, ns
 
  ! ---------------------------------------------------------

  nLev = size(ql_stoch, 3)
  nCol = size(ql_stoch, 4)
  
  !
  ! Normally, we use the overlap specified in the namelist for this module, but
  !   this can be overridden with an optional argument. 
  !
  if(present(overlap)) then
    overlapToUse = overlap
  else
    overlapToUse = defaultOverlap
  end if
  
  !
  ! Ensure that cloud fraction, water vapor, and water and ice contents are 
  ! in bounds
  !
  !   After similar code in cloud_summary3
  !
  qa_local(:,:,:) = 0.
  do k=1, nLev
       do j=1, size(qa_stoch,2)
       do i=1, size(qa_stoch,1)
         if (qa(i,j,k) > qamin .and. ql(i,j,k) > qcmin) then
           qa_local(i,j,k) = qa(i,j,k)
           ql_local(i,j,k) = ql(i,j,k)
         else
           ql_local(i,j,k) = 0.0           
         endif
         if (qa(i,j,k) > qamin .and. qi(i,j,k) > qcmin) then
           qa_local(i,j,k) = qa(i,j,k)
           qi_local(i,j,k) = qi(i,j,k)
         else
           qi_local(i,j,k) = 0.0           
         endif
       end do
       end do
  end do

  
  !
  ! Apply overlap assumption
  !
   if (overlapToUse == 2) then
      pdfRank(:, :, :, :) = genRandomOverlapSamples( qa_local(:, :, :), nCol, streams(:, :))
   else if (overlapToUse == 1) then
      pdfRank(:, :, :, :) = genMaximumOverlapSamples(qa_local(:, :, :), nCol, streams(:, :))
   else if (overlapToUse == 3) then
      pdfRank(:, :, :, :) = genMaxRanOverlapSamples( qa_local(:, :, :), nCol, streams(:, :))
   else if (overlapToUse == 4) then
      if(.not. present(pFull)) call error_mesg("cloud_generator_mod", &
                                               "Need to provide pFull when using overlap = 4", FATAL)
      !
      ! Height difference from hydrostatic equation with fixed scale height for T = 250 K
      !
      heightDifference(:, :, :nLev-1) = (log(pFull(:, :, 2:nLev)) - log(pFull(:, :, 1:nLev-1))) * &
                                         pressureScaleHeight 
      heightDifference(:, :, nLev) = heightDifference(:, :, nLev-1)
      !
      ! Overlap is weighted between max and random with parameter overlapWeighting (0 = random, 
      !    1 = max), which decreases exponentially with the separation distance. 
      !
      overlapWeighting(:, :, :) = exp( - heightDifference(:, :, :) / overlapLengthScale )
      pdfRank(:, :, :, :) = genWeightedOverlapSamples(qa_local(:, :, :),         &
                                                      overlapWeighting(:, :, :), & 
                             nCol, streams(:, :))
    else
      call error_mesg("cloud_generator_mod", "unknown overlap parameter", FATAL)
     
     endif



  
  if(.not. do_inhomogeneous_clouds .and. .not. do_pdf_clouds) then
    !
    ! The clouds are uniform, so every cloudy cell at a given level 
    !   gets the same ice and/or liquid concentration (subject to a minumum). 
    ! We're looking for in-cloud ice and water contents, so we 
    !   divide by cloud fraction
    !  
    do n=1,nCol
      do k=1,nLev
       do j=1, size(qa_stoch,2)
       do i=1, size(qa_stoch,1)
!  a "true" for the following test indicates the presence of cloudiness
         if (pdfRank(i,j,k,n) > (1.-qa_local(i,j,k))) then
           cld_thickness(i,j,k,n   ) = 1.
           qa_stoch(i,j,k,n   ) = 1. 
           ql_stoch(i,j,k,n   ) = ql_local(i,j,k)/qa_local(i,j,k)
           qi_stoch(i,j,k,n   ) = qi_local(i,j,k)/qa_local(i,j,k)
         else
           cld_thickness(i,j,k,n   ) = 0.
           qa_stoch(i,j,k,n   ) = 0. 
           ql_stoch(i,j,k,n   ) = 0.                              
           qi_stoch(i,j,k,n   ) = 0.                              
         endif
        end do
        end do
        end do
        end do
           

  end if
  
  if (do_inhomogeneous_clouds) then
    if(.not. present(pFull) .or. .not. present(temperature)) &
      call error_mesg("cloud_generator_mod",                 &
      "Need to provide pFull and temperature when using inhomogenous clouds", FATAL)
    !
    ! Assume that total water in each grid cell follows a symmetric beta distribution with  
    !   exponents p=q set in the namelist. Determine the normalized amount of condensate 
    !   (qc - qmin)/(qmax - qmin) from the incomplete beta distribution at the pdf rank. 
    !   Convert to physical units based on three equations
    !   qs_norm = (qs - qmin)/(qmax - qmin)
    !   qa = 1. - betaIncomplete(qs_norm, p, q)
    !   qc_mean/(qmax - qmin) = aThermo( p/(p+q) (1 - betaIncomplete(qs_norm, p+1, q)) - 
    !                                    qs_norm * (1 - cf) )
    !   The latter equation comes from integrating aThermo * (qtot - qsat) times a beta distribution 
    !   from qsat to qmax; see, for example, from Tompkins 2002 (Eq. 14), but including a thermodynamic 
    !   term aThermo = 1/(1 + L/cp dqs/dt) evaluated at the "frozen temperature", as per SAK. 
    !   
    where (qa_local(:, :, :) > qamin) 
      qlqcRatio(:, :, :) = ql_local(:, :, :) / (ql_local(:, :, :) + qi_local(:, :, :))
    elsewhere
      qlqcRatio(:, :, :) = 1 ! Shouldn't be used. 
    end where
    call computeThermodynamics(temperature, pFull, ql, qi, aThermo)

    qs_norm(:, :, :) = 1. ! This assignment is made so the values of qs_norm 
                          ! are always valid; in practice these should be masked out 
  ! in the regions below. 
    where(qa_local(:, :, :) < qamin)
      !
      ! Not cloudy, so deltaQ is irrelevant
      !
      qs_norm(:, :, :) = 1. 
      deltaQ(:, :, :) = 0.
    elsewhere
      !
      ! Hey, is this the right test for fully cloudy conditions? Is cloud fraction ever 1.? 
      !
      where (qa_local(:, :, :) >= 1.) 
        !
        ! In fully cloudy conditions we don't have any information about the bounds
        !   of the total water distribution. We arbitrarily set the lower bound to qsat.
        !   For a symmetric distribution the upper bound to qsat plus twice the amount of 
        !   condensate. 
        !
          qs_norm(:, :, :) = 0.
        deltaQ(:, :, :) = (2/aThermo(:, :, :)) * &
                          (ql_local(:, :, :) + qi_local(:, :, :)) / qa_local(:, :, :)
        ! qMin(:, :, :) = qsat(:, :, :)
      elsewhere 
        !
        ! Partially cloudy conditions - diagnose width of distribution from mean 
        !   condensate amount and cloud fraction. 
        !   The factor 1/2 = p/(p+q)
        !
        qs_norm(:, :, :) = beta_deviate(1. - qa_local(:, :, :), p = betaP, q = betaP)
        deltaQ(:, :, :) =                                                               &
          (ql_local(:, :, :) + qi_local(:, :, :)) /                                     &
          (aThermo(:, :, :) * ((1./2. * (1. - incomplete_beta(qs_norm(:, :, :),         &
                                                          p = betaP + 1, q = betaP))) - &
                               qs_norm(:, :, :) * qa_local(:, :, :) ))
        ! qMin(:, :, :) = qsat(:, :, :) - qs_norm(:, :, :) * deltaQ(:, :, :)
      end where 
    end where
  
    do n=1,nCol
      do k=1,nLev
       do j=1, size(qa_stoch,2)
       do i=1, size(qa_stoch,1)
!  a "true" for the following test indicates the presence of cloudiness
         if (pdfRank(i,j,k,n) > (1.-qa_local(i,j,k))) then
           cld_thickness(i,j,k,n   ) = 1.
           qa_stoch(i,j,k,n   ) = 1. 
      !
      ! Hey, do we need to account for cloud fraction here, as we do in the homogeneous case? 
      !   Also, does this seem like the right way to go from the mean to individual samples? 
      !
      qc_stoch(i, j, k, n) = aThermo(i,j,k  ) * deltaQ(i,j,k  ) * &
            (beta_deviate(pdfRank(i,j,k,n   ), p = betaP, q = betaP) - &
                   qs_norm(i,j,k) ) 
      ! 
      ! The proportion of ice and water in each sample is the same as the mean proportion.  
      !
           ql_stoch(i,j,k,n   ) = qc_stoch(i,j,k,n)* qlqcRatio(i,j,k)
           qi_stoch(i,j,k,n   ) = qc_stoch(i,j,k,n)*   &
                                     (1.-qlqcRatio(i,j,k))
         else
           cld_thickness(i,j,k,n   ) = 0.
           qa_stoch(i,j,k,n   ) = 0. 
           ql_stoch(i,j,k,n   ) = 0.                              
           qi_stoch(i,j,k,n   ) = 0.                              
         endif
        end do
        end do
        end do
        end do
           
  end if 

  if (do_pdf_clouds) then
    if(.not. present(pFull) .or. .not. present(temperature) .or. &
       .not. present(qv)    .or. .not. present(pHalf)) &

      call error_mesg("cloud_generator_mod",                 &
      "Need to provide pFull, pHalf, temperature and water vapor when using pdf clouds", FATAL)
    !
    ! Assume that total water in each grid cell follows a symmetric beta distribution with  
    !   exponents p=q set in the namelist. In contrast to the procedure for inhomogen-
    !   eous clouds, here the distribution is known (either from prognostic variance)
    !   or from diagnostic variance. In this case one can directly determine
    !   the cloud condensate from:
    !   
    !   qc = aThermo * (qt -qs) =
    !      = aThermo * deltaQ * (beta_deviate(pdfRank,p,q) - qsnorm)
    !
    !   Note that the qlqcRatio needs to be defined in the case that there is no
    !   cloud in the grid box from what enters this subroutine (qa_local) yet there
    !   is condensate diagnosed in this routine. The formula that weights the
    !   saturation vapor pressure as a function of temperature is used.
    
       
    where (qa_local(:, :, :) > qamin) 
      qlqcRatio(:, :, :) = ql_local(:, :, :) / (ql_local(:, :, :) + qi_local(:, :, :))
    elsewhere
      qlqcRatio(:, :, :) = min(1., max(0., 0.05*(temperature(:,:,:)-tfreeze+20.))) 
    end where
    call computeThermodynamics(temperature, pFull, ql, qi, aThermo, qs)

    !Compute ppm interpolation - if nsublevels > 1
    do k = 1, nLev
         delp(:,:,k) = pHalf(:,:,k+1)-pHalf(:,:,k)
    enddo     
    qta4(1,:,:,:) = max(qcmin,qv+ql_local+qi_local)
    qtqsa4(1,:,:,:) = qta4(1,:,:,:)-qs
        
    if (nsublevels.gt.1) then
        do i=1, size(qa_stoch,1)
            call ppm2m_sak(qta4(:,i,:,:),delp(i,:,:),nLev,kmap,1,&
                       size(qa_stoch,2),0,kord)
            call ppm2m_sak(qtqsa4(:,i,:,:),delp(i,:,:),nLev,kmap,1,&
                       size(qa_stoch,2),0,kord)
        enddo                
    else
        qta4(2,:,:,:) = qta4(1,:,:,:)
        qta4(3,:,:,:) = qta4(1,:,:,:)
        qta4(4,:,:,:) = 0.
        qtqsa4(2,:,:,:) = qtqsa4(1,:,:,:)
        qtqsa4(3,:,:,:) = qtqsa4(1,:,:,:)
        qtqsa4(4,:,:,:) = 0.   
    end if
    
    !loop over vertical levels       
    do k=1, nLev
         
        !initialize summing variable
        qctmp(:,:,:) = 0.         
        qatmp(:,:,:) = 0.
        
        !Create loop over sub-levels within a grid box
        do ns = 1, nsublevels
             
             !calculate normalized vertical level
             ! 0. = top of gridbox
             ! 1. = bottom of gridbox
        
             pnorm =  (real(ns) - 0.5 )/real(nsublevels)
             
             !First step is to calculating the 
             !the width of the qt distribution (deltaQ)
             deltaQ(:,:,k) = max(qcmin, &
                          qta4(2,:,:,k)+pnorm*( (qta4(3,:,:,k)- &
                          qta4(2,:,:,k)) +  qta4(4,:,:,j)*(1-pnorm) ) )
             deltaQ(:,:,k) = 2.*qthalfwidth*deltaQ(:,:,k)
             
             !From this the variable normalized saturation specific
             !humidity qs_norm is calculated.
             !
             !  qs_norm = (qs(Tl) - qtmin)/(qtmax-qtmin)
             !
             !          = 0.5  - (qtbar - qs(Tl))/deltaQ
             !
             !Note that if qs_norm > 1., the grid box is fully clear.
             !If qs_norm < 0., the grid box is fully cloudy.
             qs_norm(:,:,k) = qtqsa4(2,:,:,k)+  &
                       pnorm*( (qtqsa4(3,:,:,k)-qtqsa4(2,:,:,k)) + &
                       qtqsa4(4,:,:,k)*(1-pnorm) )
             qs_norm(:,:,k) = 0.5 - ( qs_norm(:,:,k)/deltaQ(:,:,k) )
             
             do j=1, size(qa_stoch,2)
             do i=1, size(qa_stoch,1)             
             do n=1,nCol
       
                  if (qs_norm(i,j,k).lt.1.) then
                      tmpr =           aThermo(i,j,k  ) * &
                                        deltaQ(i,j,k  ) * &
                               (beta_deviate(pdfRank(i,j,k,n   ), &
                                          p = betaP, q = betaP) - &
                                       qs_norm(i,j,k) )
                      if (tmpr.gt.qcmin) then                 
                          qctmp(i, j, n) = qctmp(i, j, n) + tmpr
                          qatmp(i, j, n) = qatmp(i, j, n) + 1.
                      end if           
                  end if
                  
             enddo !for nCol
             enddo !for i
             enddo !for j
             
        enddo !for ns (nsublevels)
         
        !produce vertical average qc
        do j=1, size(qa_stoch,2)
        do i=1, size(qa_stoch,1)
        do n=1, nCol
            
             qctmp(i,j,n) = qctmp(i,j,n) / real (nsublevels)
             qa_stoch(i,j,k,n) = qatmp(i,j,n) / real (nsublevels)
             if (qctmp(i,j,n).le.qcmin) then
                  cld_thickness(i,j,k,n) = 0.
                  qa_stoch(i,j,k,n) = 0.
                  qc_stoch(i,j,k,n) = 0.
                  ql_stoch(i,j,k,n) = 0.
                  qi_stoch(i,j,k,n) = 0.
             else
                  qc_stoch(i,j,k,n) = qctmp(i,j,n)             
                  cld_thickness(i,j,k,n) = 1.
                  !qa_stoch(i,j,k,n) = 1.          
       
                  ! 
                  ! The proportion of ice and water in each 
                  ! sample is the same as the mean proportion.  
                  !
                  ql_stoch(i,j,k,n   ) = qc_stoch(i,j,k,n)* &
                                        qlqcRatio(i,j,k)
                  qi_stoch(i,j,k,n   ) = qc_stoch(i,j,k,n)* &
                                     (1.-qlqcRatio(i,j,k))
             end if !for column containing cloud
             
        enddo !for nCol
        enddo !for i
        enddo !for j 
        
    end do !for k (vertical loop)
           
  end if  ! for do_pdf_clouds


  ! Create qn_stoch - the stochastic cloud droplet number
  ! 
  do n=1,nCol
  do k=1,nLev
  do j=1, size(qa_stoch,2)
  do i=1, size(qa_stoch,1)
    if (ql_stoch(i,j,k,n)>qcmin .and. qa_stoch(i,j,k,n)>qcmin) then
      qn_stoch(i,j,k,n) = qa_stoch(i,j,k,n)* &
                          max(0.,qn(i,j,k)/max(qcmin,qa(i,j,k)))
      !note that qa is compared to qcmin to minimize the impact of the
      !limiters on calculating the in-cloud cloud droplet number
    else
      qn_stoch(i,j,k,n) = 0.                           
    endif
  end do
  end do
  end do
  end do


!cms++
  ! Create qni_stoch - the stochastic cloud ice number
  ! 
  do n=1,nCol
  do k=1,nLev
  do j=1, size(qa_stoch,2)
  do i=1, size(qa_stoch,1)
    if (qi_stoch(i,j,k,n)>qcmin .and. qa_stoch(i,j,k,n)>qcmin) then
      qni_stoch(i,j,k,n) = qa_stoch(i,j,k,n)* &
                          max(0.,qni(i,j,k)/max(qcmin,qa(i,j,k)))
      !note that qa is compared to qcmin to minimize the impact of the
      !limiters on calculating the in-cloud cloud ice number
    else
      qni_stoch(i,j,k,n) = 0.                           
    endif
  end do
  end do
  end do
  end do
!cms--

end subroutine generate_stochastic_clouds

  ! ---------------------------------------------------------
  !  Function to return the weighting between maximum and random overlap 
  !    given the pressure difference 
  ! 
  !  Note pPlus is the pressure at a higher altitude
  !  (i.e. pPlus < pMinus)
  ! ---------------------------------------------------------
function compute_overlap_weighting(qaPlus, qaMinus, pPlus, pMinus) result(weighting)
  real, dimension(:, :), intent( in) :: qaPlus, qaMinus, pPlus, pMinus
  real, dimension(size(pPlus, 1), &
                  size(pPlus, 2))    :: weighting
        
  select case(defaultOverlap)
    case(1) ! Maximum overlap
      weighting(:, :) = 1.
    case(2) ! Random overlap
      weighting(:, :) = 0.
    case(3) ! Maximum-random
      where(qaPlus(:, :) > qamin) 
        weighting(:, :) = 1.
      elsewhere
        weighting(:, :) = 0.
      end where
    case(4)
      !
      ! Overlap is weighted between max and random with parameter overlapWeighting (0 = random, 
      !    1 = max), which decreases exponentially with the separation distance. 
      !
      weighting(:, :) = exp(-abs(log(pMinus(:, :)) - log(pPlus(:, :))) * &
                              pressureScaleHeight / overlapLengthScale)
    case default 
      call error_mesg("cloud_generator_mod", "unknown overlap parameter", FATAL)
    end select
      
end function compute_overlap_weighting



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!  These generate cloud samples according to specific overlap rules.  
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine computeThermodynamics(temperature, pressure, ql, qi, aThermo, qs)
    real, dimension(:, :, :), intent( in) :: temperature, pressure, ql, qi
    real, dimension(:, :, :), intent(out) :: aThermo
    real, dimension(:, :, :), optional, intent(out) :: qs
        
    real, parameter :: d608 = (rvgas - rdgas)/rdgas, &
                       d622 = rdgas/rvgas,           &
                       d378 = 1. - d622
    
    ! Compute saturation mixing ratio and gamma = 1/(1 + L/cp dqs/dT) evaluated
    !   at the ice water temperature
    ! Taken from strat_cloud_mod
    
    ! Local variables
    real, dimension(size(temperature, 1), &
                    size(temperature, 2), &
                    size(temperature, 3)) :: Tl, L, dqsdT,qsloc
    
    !
    ! Ice water temperature - ql and qi are grid cell means
    !
    Tl(:, :, :) =  temperature(:, :, :) -       &
                   (hlv/cp_air) * ql(:, :, :) - &
                   (hls/cp_air) * qi(:, :, :)
  
    !calculate qs and dqsdT
    call compute_qs (Tl, pressure,  qsloc, dqsdT = dqsdT)
    if (present(qs)) then
       qs = qsloc
    endif
 
    ! Latent heat of phase change, varying from that of water to ice with temperature.
    ! 
    L(:, :, :) = (min(1., max(0., 0.05*(temperature(:, :, :) - tfreeze + 20.)))*hlv + &
                  min(1., max(0., 0.05*(tfreeze - temperature(:, :, :)      )))*hls)
    aThermo(:, :, :) = 1./ (1. + L(:, :, :)/cp_air * dqsdT(:, :, :)) 
    
  end subroutine computeThermodynamics
  ! ---------------------------------------------------------
  ! Random overlap - the value is chosen randomly from the distribution at 
  !   every height. 
  ! ---------------------------------------------------------
  function genRandomOverlapSamples_1D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,    dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    ! -------------------
    call getRandomNumbers(stream, randomValues)
  end  function genRandomOverlapSamples_1D
  ! ---------------------------------------------------------
  function genRandomOverlapSamples_3D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :),    intent(in   ) :: cloudFraction
    integer,                        intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),       intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)                  :: randomValues
    ! Local variables
    integer :: i, j
    ! -------------------
    do j = 1, size(cloudFraction, 2)
      do i = 1, size(cloudFraction, 1)
        call getRandomNumbers(stream(i, j), randomValues(i, j, :, :))
      end do
    end do
  end  function genRandomOverlapSamples_3D
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  ! Maximum overlap - the position in the PDF is the same 
  !   at every height in a given column (though it varies from 
  !   column to column). 
  ! ---------------------------------------------------------

  function genMaximumOverlapSamples_1D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,     dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    
    
    ! -------------------
    call getRandomNumbers(stream, randomValues(1, :))
    randomValues(:, :)  = spread(randomValues(1, :), &
                                 dim = 1, nCopies = size(cloudFraction))
  end  function genMaximumOverlapSamples_1D
  ! ---------------------------------------------------------
  function genMaximumOverlapSamples_3D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :),    intent(in   ) :: cloudFraction
    integer,                        intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),       intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)                  :: randomValues
    ! -------------------
    ! Local variables
    integer :: i, j, nX, nY, nLev
    ! -------------------
    nX   = size(cloudFraction, 1); nY  = size(cloudFraction, 2)
    nLev = size(cloudFraction, 3)
    
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues(i, j, 1, :))
      end do
    end do 
    randomValues(:, :, :, :)  = spread(randomValues(:, :, 1, :), dim = 3, nCopies = nLev)

  end  function genMaximumOverlapSamples_3D
  ! ---------------------------------------------------------
  
  ! ---------------------------------------------------------
  ! Meximum-random overlap. 
  ! Within each column, the value in the top layer is chosen 
  !   at random. We then walk down one layer at a time. If the layer above is cloudy
  !   we use the same random deviate in this layer; otherwise
  !   we choose a new value. 
  ! ---------------------------------------------------------
  
  function genMaxRanOverlapSamples_1D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,    dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    ! Local variables
    integer                                 :: level
    
    ! -------------------
    call getRandomNumbers(stream, randomValues)
    do level = 2, size(cloudFraction)
      where(randomValues(:, level - 1) > 1. - cloudFraction(level - 1))
        randomValues(level, :) = randomValues(level - 1, :)
      elsewhere
        randomValues(level, :) = randomValues(level, :) * (1. - cloudFraction(level - 1))
      end where
    end do
  end  function genMaxRanOverlapSamples_1D
  ! ---------------------------------------------------------
  function genMaxRanOverlapSamples_3D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :),    intent(in   ) :: cloudFraction
    integer,                        intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),       intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)                  :: randomValues
    ! -------------------
    ! Local variables
    integer :: i, j, level, nX, nY, nLev
    ! -------------------
    nX   = size(cloudFraction, 1); nY  = size(cloudFraction, 2)
    nLev = size(cloudFraction, 3)
    
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues(i, j, :, :))
      end do
    end do 
              
    do level = 2, nLev
      where(randomValues(:, :, level - 1, :) > &
            spread(1. - cloudFraction(:, :, level - 1), dim = 3, nCopies = nSamples))
        randomValues(:, :, level, :) = randomValues(:, :, level - 1, :)
      elsewhere
        randomValues(:, :, level, :) = randomValues(:, :, level, :) * &
                                       spread(1. - cloudFraction(:, :, level - 1), &
                                              dim = 3, nCopies = nSamples)
      end where
    end do
  end  function genMaxRanOverlapSamples_3D
  ! ---------------------------------------------------------
  
  ! Neighbor to neighbor rank correlation, which gives exponential dependence 
  !   of rank correlation if the correlation is fixed. The correlation coefficient 
  !   is the array alpha. 
  ! Two streams of random numbers are generated. The first corresponds to the postion in 
  !   the PDF, and the second is used to enforce the correlation . If the value of 
  !   the second stream at one level in one column is less than alpha at that level, 
  !   the same relative position in the PDF is chosen in the lower layer as the upper.  
  ! ---------------------------------------------------------
  function genWeightedOverlapSamples_1D(cloudFraction, alpha, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction, alpha
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,    dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    
    ! Local variables
    real, dimension(size(cloudFraction), nSamples) :: randomValues2
    integer                                        :: level
    
    ! -------------------
    call getRandomNumbers(stream, randomValues)
    call getRandomNumbers(stream, randomValues2)
    
    do level = 1, size(cloudFraction) - 1 
      where(randomValues2(level + 1, :) < alpha(level)) &
        randomValues(level + 1, :) = randomValues(level, :)
    end do   
  end  function genWeightedOverlapSamples_1D
  ! ---------------------------------------------------------
  function genWeightedOverlapSamples_3D(cloudFraction, alpha, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :), intent(in   ) :: cloudFraction, alpha
    integer,                     intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),    intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)               :: randomValues
     
    ! Local variables
    real, dimension(size(cloudFraction, 1), &
                    size(cloudFraction, 2), &
                    size(cloudFraction, 3), &
                                  nSamples) :: randomValues2
    integer                                 :: i, j, nX, nY, nLev, level
    
    ! -------------------
    nX   = size(cloudFraction, 1); nY  = size(cloudFraction, 2)
    nLev = size(cloudFraction, 3)
    
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues (i, j, :, :))
      end do
    end do 
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues2(i, j, :, :))
      end do
    end do 
     
    do level = 1, nLev - 1
      where(randomValues2(:, :, level + 1, :) < spread(alpha(:, :, level),           &
                                                       dim = 3, nCopies = nSamples)) &
        randomValues(:, :, level + 1, :) = randomValues(:, :, level, :)
    end do
  end  function genWeightedOverlapSamples_3D
  ! ---------------------------------------------------------

  subroutine cloud_generator_end       
  !----------------------------------------------------------------------
  !    cloud_generator_end is the destructor for cloud_generator_mod.
  !----------------------------------------------------------------------
          
  !---------------------------------------------------------------------
  !    be sure module has been initialized.
  !---------------------------------------------------------------------
        if (.not. module_is_initialized ) then
          call error_mesg ('cloud_generator_mod',   &
               'module has not been initialized', FATAL )
        endif
        
        if(do_inhomogeneous_clouds .or. do_pdf_clouds) &
             call beta_dist_end
  !---------------------------------------------------------------------
  !    mark the module as not initialized.
  !---------------------------------------------------------------------
        module_is_initialized = .false.
  end subroutine cloud_generator_end
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !
  !  Function to report if the cloud generator is being used. 
  !
  function do_cloud_generator()
    logical :: do_cloud_generator
    
    do_cloud_generator = cloud_generator_on
  end function do_cloud_generator
  !--------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m_sak --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m_sak(a4, delp, km, kmap, i1, i2, iv, kord)

 implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kmap    ! partial remap to start
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real, intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real   dc(i1:i2,km)
      real   h2(i1:i2,km)
      real delq(i1:i2,km)
      real  df2(i1:i2,km)
      real   d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt
      integer it
      real fac
      real a1, a2, c1, c2, c3, d1, d2
      real qmax, qmin, cmax, cmin
      real qm, dq, tmp
      real qmp, pmp
      real lac

      km1 = km - 1
       it = i2 - i1 + 1

      do k=max(2,kmap-2),km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo
 
      do k=max(2,kmap-2),km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=max(3,kmap), km1
      do i=i1,i2
        c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                  ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      if(km>8 .and. kord>3) call steepz_sak(i1, i2, km, kmap, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      if ( kmap <= 2 ) then
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
         dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
         cmax = max(a4(1,i,1), a4(1,i,2))
         cmin = min(a4(1,i,1), a4(1,i,2))
         a4(2,i,2) = max(cmin,a4(2,i,2))
         a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

      do k=max(1,kmap),km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo

! Enforce monotonicity of the "slope" within the top layer
      if ( kmap <= 2 ) then
      do i=i1,i2
         if ( a4(2,i,1) * a4(1,i,1) <= 0. ) then 
              a4(2,i,1) = 0.
                dc(i,1) = a4(1,i,1)
         endif
         if ( dc(i,1) * (a4(2,i,2) - a4(1,i,1)) <= 0. ) then
! Setting DC==0 will force piecewise constant distribution after
! calling kmppm_sak
              dc(i,1) = 0.
         endif
      enddo
      endif

! Enforce constraint on the "slope" at the surface

      do i=i1,i2
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) then
!            a4(3,i,km) = 0.
!              dc(i,km) =  -a4(1,i,km)
               dc(i,km) = 0.
         endif
         if( dc(i,km) * (a4(1,i,km) - a4(2,i,km)) <= 0. ) then
             dc(i,km) = 0.
         endif
      enddo
 
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      if ( kmap <= 2 ) then
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
            call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
      endif

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=max(2,kmap-1), km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=max(3,kmap), km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=max(3,kmap), km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
!EOC
 end subroutine ppm2m_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm_sak --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm_sak(dm, a4, itot, lmt)

 implicit none

! !INPUT PARAMETERS:
      real, intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real, parameter:: r12 = 1./12.
      real qmp
      real da1, da2, a6da
      real fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2003)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

!EOC
 end subroutine kmppm_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz_sak --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz_sak(i1, i2, km, kmap, a4, df2, dm, dq, dp, d4)

   implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: kmap                 ! 
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real, intent(in) ::  dp(i1:i2,km)       ! grid size
      real, intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real, intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real, intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real, intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real alfa(i1:i2,km)
      real    f(i1:i2,km)
      real  rat(i1:i2,km)
      real  dg2

! Compute ratio of dq/dp
      do k=max(2,kmap-1),km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=max(2,kmap-1),km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=max(3,kmap),km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=max(4,kmap+1),km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

!EOC
 end subroutine steepz_sak
!----------------------------------------------------------------------- 

end module cloud_generator_mod

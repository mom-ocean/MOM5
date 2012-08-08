module beta_dist_mod
  use fms_mod,only: error_mesg, FATAL, WARNING
  
  implicit none
  private 
  
  ! Provide values of the beta distribtion as a function of the CDF (the incomplete beta
  !   function). Returns a value as a function of two beta distribution parameters p, q 
  !   (here they must be integers) and the value x of the CDF between 0 and 1. 
  
  ! In this version we build tables using the NAG library function nag_beta_deviate, then 
  !   look up values from a table. The table can be built at run time or read from 
  !   a file (this version uses netcdf format). 
  
  ! betaDeviateTable is a 3D table with dimensions
  !   x, p, q. The range of P and Q are specified when the tables are built. 
  !   The arrays bounds are from 0 to nSteps + 1, just in case we draw exactly 0 or 1. 
  !
  character(len=128)  :: version =  '$Id: betaDistribution.F90,v 13.0 2006/03/28 21:07:32 fms Exp $'
  character(len=128)  :: tagname =  '$Name: siena_201207 $'
  
  logical         :: module_is_initialized = .false.
  

  interface beta_deviate
    module procedure betaDeviate_s,  betaDeviate_1D, betaDeviate_2D, &
                     betaDeviate_3D, betaDeviate_4D
  end interface ! betaDeviate

  interface incomplete_beta
    module procedure incompleteBeta_s,  incompleteBeta_1D, incompleteBeta_2D, &
                     incompleteBeta_3D, incompleteBeta_4D
  end interface ! incompleteBeta

  public :: beta_dist_init, beta_deviate, incomplete_beta, beta_dist_end
contains
 ! ---------------------------------------------------------
  subroutine test_beta
  
    integer :: i
    real    :: x, inc_x, inv_inc_x, inv_x, inc_inv_x
    
  end subroutine test_beta
  ! ---------------------------------------------------------
  subroutine beta_dist_init
    ! Initialize the tables containing the incomplete beta function
    !   and its inverse (beta deviate). 
    !   If the table parameters are supplied we use the NAG libraries to 
    !   compute a new table and write it to a file; if just
    !   the file name is supplied we read the table from the file. 
    !
    
    call error_mesg('beta_dist_mod', &
      'This module is not supported as part of the public release', FATAL)
    
  end subroutine beta_dist_init
!---------------------------------------------------------------------
  subroutine beta_dist_end

    !---------------------------------------------------------------------
    !    be sure module has been initialized.
    !---------------------------------------------------------------------
    call error_mesg('beta_dist_mod', &
      'This module is not supported as part of the public release', WARNING)
    
  end subroutine beta_dist_end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     SEMI-PRIVATE PROCEDURES 
!            Not accessed directly but through generic interface
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ---------------------------------------------------------
!  Functions to look up the beta deviate (inverse incomplete beta distribution) 
!    from a table
!    Overloaded, to allow for input arguments from 0 to 4 dimensions
!    It might be more efficient to loop over dimensions higher than 1 to 
!    avoid using reshape.
! ---------------------------------------------------------
  function betaDeviate_s(x, p, q) result (betaDeviate)
    real,                  intent( in) :: x
    integer,               intent( in) :: p, q
    real                               :: betaDeviate
    
    call error_mesg('betaDeviate', &
      'This module is not supported as part of the public release', FATAL)
    
    betaDeviate = 0.0
  end function betaDeviate_s
! ---------------------------------------------------------
  function betaDeviate_1D(x, p, q) result (betaDeviate)
    real,    dimension(:),    intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x))           :: betaDeviate
    
    call error_mesg('betaDeviate', &
      'This module is not supported as part of the public release', FATAL)
    
    betaDeviate(:) = 0.0

  end function betaDeviate_1D
! ---------------------------------------------------------
  function betaDeviate_2D(x, p, q) result (betaDeviate)
    real,    dimension(:, :), intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2))        :: betaDeviate
    
    call error_mesg('betaDeviate', &
      'This module is not supported as part of the public release', FATAL)
    
    betaDeviate(:, :) = 0.0
    
  end function betaDeviate_2D
! ---------------------------------------------------------
  function betaDeviate_3D(x, p, q) result (betaDeviate)
    real,    dimension(:, :, :), intent( in) :: x
    integer,                     intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2), &
                       size(x, 3))           :: betaDeviate
    
    call error_mesg('betaDeviate', &
      'This module is not supported as part of the public release', FATAL)
    
    betaDeviate(:, :, :) = 0.0

  end function betaDeviate_3D
! ---------------------------------------------------------
  function betaDeviate_4D(x, p, q) result (betaDeviate)
    real,    dimension(:, :, :, :), intent( in) :: x
    integer,                        intent( in) :: p, q
    real,    dimension(size(x, 1), size(x, 2), &
                       size(x, 3), size(x, 4))  :: betaDeviate
    
    call error_mesg('betaDeviate', &
      'This module is not supported as part of the public release', FATAL)
    
    betaDeviate(:, :, :, :) = 0.0

  end function betaDeviate_4D
! ---------------------------------------------------------
! ---------------------------------------------------------
!  Functions to look up the incomplete beta function from a table. 
!    Overloaded, to allow for input arguments from 0 to 4 dimensions
!    It might be more efficient to loop over dimensions higher than 1 to 
!    avoid using reshape.
! ---------------------------------------------------------
  function incompleteBeta_s(x, p, q) result (incompleteBeta)
    real,                  intent( in) :: x
    integer,               intent( in) :: p, q
    real                               :: incompleteBeta
    
    call error_mesg('incompleteBeta', &
      'This module is not supported as part of the public release', FATAL)
    
    incompleteBeta = 0.0

  end function incompleteBeta_s
! ---------------------------------------------------------
  function incompleteBeta_1D(x, p, q) result (incompleteBeta)
    real,    dimension(:),    intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x))           :: incompleteBeta
    
    call error_mesg('incompleteBeta', &
      'This module is not supported as part of the public release', FATAL)
    
    incompleteBeta(:) = 0.0

  end function incompleteBeta_1D
! ---------------------------------------------------------
  function incompleteBeta_2D(x, p, q) result (incompleteBeta)
    real,    dimension(:, :), intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2))        :: incompleteBeta
    
    call error_mesg('incompleteBeta', &
      'This module is not supported as part of the public release', FATAL)
    
    incompleteBeta(:, :) = 0.0

  end function incompleteBeta_2D
! ---------------------------------------------------------
  function incompleteBeta_3D(x, p, q) result (incompleteBeta)
    real,    dimension(:, :, :), intent( in) :: x
    integer,                     intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2), &
                       size(x, 3))           :: incompleteBeta
    
    call error_mesg('incompleteBeta', &
      'This module is not supported as part of the public release', FATAL)
    
    incompleteBeta(:, :, :) = 0.0
  end function incompleteBeta_3D
! ---------------------------------------------------------
  function incompleteBeta_4D(x, p, q) result (incompleteBeta)
    real,    dimension(:, :, :, :), intent( in) :: x
    integer,                        intent( in) :: p, q
    real,    dimension(size(x, 1), size(x, 2), &
                       size(x, 3), size(x, 4))  :: incompleteBeta
    
    call error_mesg('incompleteBeta', &
      'This module is not supported as part of the public release', FATAL)
    
    incompleteBeta(:, :, :, :) = 0.0

  end function incompleteBeta_4D
! ---------------------------------------------------------
end module beta_dist_mod


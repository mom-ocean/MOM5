module beta_dist_mod
  use fms_mod,only: stdlog, write_version_number, &
                    error_mesg, FATAL
  use mpp_mod,only: get_unit                  
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
  character(len=128)  :: version =  '$Id: betaDistribution.F90,v 16.0 2008/07/30 22:06:18 fms Exp $'
  character(len=128)  :: tagname =  '$Name: tikal $'
  
  logical         :: module_is_initialized = .false.
  
  real, parameter :: failureValue = -1. 
  real, parameter :: Xmin = 0., Xmax = 1.

  integer         :: Pmax, Qmax, numXSteps
  real, dimension(:, :, :), allocatable &
                  :: betaDeviateTable, incompleteBetaTable
  character(len = 32) &
                  :: fileName = "INPUT/BetaDistributionTable.txt" 
  
  interface interpolateFromTable
    module procedure interpolateFromTable_s,  interpolateFromTable_1D, interpolateFromTable_2D, &
                     interpolateFromTable_3D, interpolateFromTable_4D
  end interface ! interpolateFromTable

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
    
    print *, "TESTING BETA" 
    print *, "x    inc(x)    inv(inc(x))   inv(x)   inc(inv(x))"
    do i = 0, numXSteps
      x = real(i)/real(numXSteps) 
      inc_x     = incomplete_beta(    x, 5, 5)
      inv_inc_x = beta_deviate(   inc_x, 5, 5)
      inv_x     = beta_deviate(       x, 5, 5)
      inc_inv_x = incomplete_beta(inv_x, 5, 5)
      write(*, "(5(f10.7, 1x))") x, inc_x, inv_inc_x, inv_x, inc_inv_x
    end do
  end subroutine test_beta
  ! ---------------------------------------------------------
  subroutine beta_dist_init
    ! Initialize the tables containing the incomplete beta function
    !   and its inverse (beta deviate). 
    !   If the table parameters are supplied we use the NAG libraries to 
    !   compute a new table and write it to a file; if just
    !   the file name is supplied we read the table from the file. 
    !
    
    !---------------------------------------------------------------------
    !    if routine has already been executed, exit.
    !---------------------------------------------------------------------
      if (module_is_initialized) return
    
    !---------------------------------------------------------------------
    !    mark the module as initialized if we're able to read the tables
    !---------------------------------------------------------------------
    call write_version_number (version, tagname)
    module_is_initialized = readFromFile(fileName)
    
  end subroutine beta_dist_init
!---------------------------------------------------------------------
  subroutine beta_dist_end

    !---------------------------------------------------------------------
    !    be sure module has been initialized.
    !---------------------------------------------------------------------
    if (.not. module_is_initialized ) return
    
    if(allocated(betaDeviateTable))    deallocate(betaDeviateTable)
    if(allocated(incompleteBetaTable)) deallocate(incompleteBetaTable)
    module_is_initialized = .false.
    
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
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax, x < Xmin, x > Xmax /) )) then
      betaDeviate = FailureValue
    else
      betaDeviate = interpolateFromTable(x, p, q, betaDeviateTable)
    end if
  end function betaDeviate_s
! ---------------------------------------------------------
  function betaDeviate_1D(x, p, q) result (betaDeviate)
    real,    dimension(:),    intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x))           :: betaDeviate
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      betaDeviate(:) = FailureValue
    else
      betaDeviate(:) = interpolateFromTable(x, p, q, betaDeviateTable)
    end if
  end function betaDeviate_1D
! ---------------------------------------------------------
  function betaDeviate_2D(x, p, q) result (betaDeviate)
    real,    dimension(:, :), intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2))        :: betaDeviate
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      betaDeviate(:, :) = FailureValue
    else
      betaDeviate(:, :) = interpolateFromTable(x, p, q, betaDeviateTable)
    end if
  end function betaDeviate_2D
! ---------------------------------------------------------
  function betaDeviate_3D(x, p, q) result (betaDeviate)
    real,    dimension(:, :, :), intent( in) :: x
    integer,                     intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2), &
                       size(x, 3))           :: betaDeviate
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      betaDeviate(:, :, :) = FailureValue
    else
      betaDeviate(:, :, :) = interpolateFromTable(x, p, q, betaDeviateTable)
    end if
  end function betaDeviate_3D
! ---------------------------------------------------------
  function betaDeviate_4D(x, p, q) result (betaDeviate)
    real,    dimension(:, :, :, :), intent( in) :: x
    integer,                        intent( in) :: p, q
    real,    dimension(size(x, 1), size(x, 2), &
                       size(x, 3), size(x, 4))  :: betaDeviate
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      betaDeviate(:, :, :, :) = FailureValue
    else
      betaDeviate(:, :, :, :) = interpolateFromTable(x, p, q, betaDeviateTable)
    end if
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
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax, x < Xmin, x > Xmax /) )) then
      incompleteBeta = FailureValue
    else
      incompleteBeta = interpolateFromTable(x, p, q, incompleteBetaTable)
    end if
  end function incompleteBeta_s
! ---------------------------------------------------------
  function incompleteBeta_1D(x, p, q) result (incompleteBeta)
    real,    dimension(:),    intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x))           :: incompleteBeta
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      incompleteBeta(:) = FailureValue
    else
      incompleteBeta(:) = interpolateFromTable(x, p, q, incompleteBetaTable)
    end if
  end function incompleteBeta_1D
! ---------------------------------------------------------
  function incompleteBeta_2D(x, p, q) result (incompleteBeta)
    real,    dimension(:, :), intent( in) :: x
    integer,                  intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2))        :: incompleteBeta
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      incompleteBeta(:, :) = FailureValue
    else
      incompleteBeta(:, :) = interpolateFromTable(x, p, q, incompleteBetaTable)
    end if
  end function incompleteBeta_2D
! ---------------------------------------------------------
  function incompleteBeta_3D(x, p, q) result (incompleteBeta)
    real,    dimension(:, :, :), intent( in) :: x
    integer,                     intent( in) :: p, q
    real,    dimension(size(x, 1), &
                       size(x, 2), &
                       size(x, 3))           :: incompleteBeta
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      incompleteBeta(:, :, :) = FailureValue
    else
      incompleteBeta(:, :, :) = interpolateFromTable(x, p, q, incompleteBetaTable)
    end if
  end function incompleteBeta_3D
! ---------------------------------------------------------
  function incompleteBeta_4D(x, p, q) result (incompleteBeta)
    real,    dimension(:, :, :, :), intent( in) :: x
    integer,                        intent( in) :: p, q
    real,    dimension(size(x, 1), size(x, 2), &
                       size(x, 3), size(x, 4))  :: incompleteBeta
    
    if (.not. module_is_initialized ) then
      call error_mesg ('beta_dist_mod', 'module has not been initialized', FATAL )
    endif
    
    if(any( (/ p < 1, p > pMax, q < 1, q > qMax /) )) then
      incompleteBeta(:, :, :, :) = FailureValue
    else
      incompleteBeta(:, :, :, :) = interpolateFromTable(x, p, q, incompleteBetaTable)
    end if
  end function incompleteBeta_4D
! ---------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE PROCEDURES 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ---------------------------------------------------------
  function interpolateFromTable_s(x, p, q, table) result (values)
    real,                      intent( in) :: x
    integer,                   intent( in) :: p, q
    real, dimension(0:, :, :), intent( in) :: table
    real                                  :: values
    
    ! Local variables
    integer :: xIndex
    real    :: xWeight
    
    ! Check to see if we're at an endpoint of the distribution 
    if      (abs(x - xMax) < spacing(Xmax)) then
      values = 1.
    else if (abs(x - XMin) < spacing(Xmin))  then
      values = 0.
    else
      !
      ! Linear interpolation in x
      !
      xIndex = int( x * numXSteps)
      xWeight = numXSteps *  x - xIndex
      values = (1. - xWeight) * table(xIndex    , p, q) + & 
               (     xWeight) * table(xIndex + 1, p, q)
    end if
  end function interpolateFromTable_s
  ! ---------------------------------------------------------
  function interpolateFromTable_1D(x, p, q, table) result(values)
    real,    dimension(:),        intent( in) :: x
    integer,                      intent( in) :: p, q
    real,    dimension(0:, :, :), intent( in) :: table
    real,    dimension(size(x))           :: values
    
    ! Local variables
    integer, dimension(size(x)) :: xIndex
    real,    dimension(size(x)) :: xWeight
    
    ! Check for parameters out of bounds, and be sure the table has been initialized. 
    ! 
    where(x(:) < Xmin .or. x(:) > Xmax) 
      values(:) = FailureValue
    else where
      !
      ! Linear interpolation in x
      !
      xIndex(:)  = int( x(:) * numXSteps)
      xWeight(:) = numXSteps *  x(:) - xIndex(:)
      values(:)  = (1. - xWeight(:)) * table(xIndex(:)    , p, q) + & 
                   (     xWeight(:)) * table(xIndex(:) + 1, p, q)
    end where

     ! Check to see if we're at an endpoint of the distribution 
    where (abs(x(:) - xMax) < spacing(Xmax)) values(:) = 1.
    where (abs(x(:) - XMin) < spacing(Xmin)) values(:) = 0.
  end function interpolateFromTable_1D
  ! ---------------------------------------------------------
  function interpolateFromTable_2D(x, p, q, table) result(values)
    real,    dimension(:, :),     intent( in) :: x
    integer,                      intent( in) :: p, q
    real,    dimension(0:, :, :), intent( in) :: table
    real,    dimension(size(x, 1), &
                       size(x, 2))           :: values
    ! Local variables 
    integer :: i
 
    do i = 1, size(x, 2)
      values(:, i) = interpolateFromTable(x(:, i), p, q, table)
    end do 
  end function interpolateFromTable_2D
  ! ---------------------------------------------------------
  function interpolateFromTable_3D(x, p, q, table) result(values)
    real,    dimension(:, :, :),  intent( in) :: x
    integer,                      intent( in) :: p, q
    real,    dimension(0:, :, :), intent( in) :: table
    real,    dimension(size(x, 1), &
                       size(x, 2), &  
                       size(x, 3))           :: values
    ! Local variables 
    integer :: i
 
    do i = 1, size(x, 3)
      values(:, :, i) = interpolateFromTable(x(:, :, i), p, q, table)
    end do 
  end function interpolateFromTable_3D
  ! ---------------------------------------------------------
  function interpolateFromTable_4D(x, p, q, table) result(values)
    real,    dimension(:, :, :, :),  intent( in) :: x
    integer,                         intent( in) :: p, q
    real,    dimension(0:, :, :),    intent( in) :: table
    real,    dimension(size(x, 1), &
                       size(x, 2), &
                       size(x, 3), &
                       size(x, 4))              :: values
    ! Local variables 
    integer :: i

    do i = 1, size(x, 4)
      values(:, :, :, i) = interpolateFromTable(x(:, :, :, i), p, q, table)
    end do 
  end function interpolateFromTable_4D
  ! ---------------------------------------------------------
  ! Reading and writing procedures
  ! ---------------------------------------------------------
  function readFromFile(fileName)
    character(len = *), intent( in) :: fileName
    logical                         :: readFromFile

    ! Local variables
!    integer, parameter :: unit = 909
    integer :: unit
    
    unit = get_unit()
    open(unit = unit, file = trim(fileName), status = 'old')
    read(unit, '(3(i5, 1x))') Pmax, Qmax, numXSteps
    allocate(   betaDeviateTable(0:numXSteps + 1, Pmax, Qmax), &
             incompleteBetaTable(0:numXSteps + 1, Pmax, Qmax))
    read(unit, '(8(f10.8, 1x))') betaDeviateTable
    read(unit, '(8(f10.8, 1x))') incompleteBetaTable
    close(unit)
 
!    print *, "Reading beta distribution tables..."
!    write(*, '(8(f10.8, 1x))') betaDeviateTable(:8, 5, 5)
!    write(*, '(8(f10.8, 1x))') incompleteBetaTable(:8, 5, 5)
   
    readFromFile = .true.
  end function readFromFile 
  
!   function readFromFile(fileName) 
!     use netcdf
!     character(len = *), intent( in) :: fileName
!     logical                         :: readFromFile
!     
!     ! Local variables - all related to netcdf
!     integer, dimension(16) :: status
!     integer                :: fileId, ncDimId, ncVarId
!    
!     status( :) = nf90_NoErr
!     status( 1) = nf90_open(trim(fileName), nf90_NoWrite, fileID)
!     status( 2) = nf90_inq_dimid(fileId, "x", ncDimId)
!     status( 3) = nf90_Inquire_Dimension(fileId, ncDimId, len = numXSteps)
!     
!     if(allocated(betaDeviateTable)) deallocate(betaDeviateTable)
!     allocate(betaDeviateTable(0:numXSteps-1, Pmax, Qmax))
!     status( 4) = nf90_inq_varId(fileId, "betaDeviateTable", ncVarId)
!     status( 5) = nf90_get_var(fileId, ncVarId, betaDeviateTable)
!     
!     if(allocated(incompleteBetaTable)) deallocate(incompleteBetaTable)
!     allocate(incompleteBetaTable(0:numXSteps-1, Pmax, Qmax))
!     status( 7) = nf90_inq_varId(fileId, "incompleteBetaTable", ncVarId)
!     status( 8) = nf90_get_var(fileId, ncVarId, betaDeviateTable)
!  
!     status( 9) = nf90_close(fileId)
!     
!     !
!     ! Make the definitions of numXSteps, consistent with what's in routine 
!     !   initializeIncBetaTables, which has a zero at one end an one extra element in each 
!     !   direction, just in case we hit the last element. 
!     ! 
!     numXSteps = numXSteps - 2
!    
!     readFromFile = all(status(:) == nf90_NoErr)
!   end function readFromFile
end module beta_dist_mod


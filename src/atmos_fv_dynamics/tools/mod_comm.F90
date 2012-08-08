!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mod_comm --- SPMD parallel decompostion/communication module
module mod_comm
!
! !USES:
#if defined ( SPMD )

#if defined( TIMING )
      use timingModule, only : timing_on, timing_off
#endif
      use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4,  &
                               i8 => shr_kind_i8, i4 => shr_kind_i4
      implicit none

#if defined ( USE_SHMEM )

#if defined(USE_VT)
#include "mpif.h"
#endif

#if defined(IRIX64)
#include "mpp/shmem.fh"
#else
#include "shmem.fh"
#endif
#else
#include "mpif.h"
#if defined(LINUX) && (USE_MPI_IO)
#include "mpiof.h"
#endif
#endif

! !PUBLIC MEMBER FUNCTIONS:

      public mp_init, mp_exit, y_decomp, set_decomp,               &
             mp_send4d_ns, mp_recv4d_ns, mp_send2_ns, mp_recv2_ns, &
             mp_send3d, mp_recv3d, mp_send3d_2, mp_recv3d_2,       &
             mp_reduce_sum, mp_reduce_max, mp_minmax, mp_sum1d, mp_add1d, &
             mp_bcst_real, mp_bcst_int, mp_bcst_r2d, mp_gath_r2d, mp_gath_3d, &
             mp_scat3d_int, mp_scat3d, mp_gath4_r2d, &
#if defined(USE_MPI_IO)
             mp_scatter3d, mp_gather3d, mp_wrt3d,                     &
             mp_file_open, mp_file_write_at, mp_file_close
#else
             mp_scatter3d, mp_gather3d, mp_wrt3d
#endif

! !PUBLIC DATA MEMBERS:
      integer, SAVE:: gid                         ! PE id
      integer(i4), SAVE, private :: numpro        ! Permanent No. of PEs
      integer(i4), SAVE, private :: numcpu        ! No. of threads

!------------------------------------------------------------------------------
!  Local parameters for use with MPI-2 and MPI-1
!------------------------------------------------------------------------------
      integer, parameter:: nbuf = 2               ! Max No. of sends per call
      integer, parameter:: nghost = 3             ! Max No. of ghost indices
      integer, parameter:: max_call = 2           ! Max No. of back-to-back...
                                                  ! ...mp_send calls

#if defined(MPI2) || defined(USE_SHMEM)
! Statically define the model dimension
#     define PLON        144
#     define PLAT         91
#     define PLEV         48
      integer, parameter:: idimsize = PLON*nghost*(PLEV+1)
#else
! Dynamic memory allocation with MPI1
      integer, save:: PLON
      integer, save:: PLAT
      integer, save:: PLEV
      integer, save:: idimsize
#endif

#if defined(USE_SHMEM)
      integer, parameter :: mp_r4 = r4
      integer, parameter :: mp_r8 = r8
      integer, parameter :: mp_i4 = i4
      integer, SAVE::  reduce_sync(SHMEM_REDUCE_SYNC_SIZE)
      integer, SAVE:: collect_sync(SHMEM_COLLECT_SYNC_SIZE)
      integer, SAVE::   bcast_sync(SHMEM_BCAST_SYNC_SIZE)
#else
      integer, parameter :: mp_r4 = MPI_REAL
      integer, parameter :: mp_r8 = MPI_DOUBLE_PRECISION
      integer, parameter :: mp_i4 = MPI_INTEGER
#endif

!------------------------------------------------------------------------------
!  Local variables for use with MPI-2, MPI-1, and SHMEM
!------------------------------------------------------------------------------
      integer, SAVE:: gsize                       ! No. of PEs
      integer, SAVE:: nowpro                      ! Temp. PE id
      integer, SAVE:: np_loop                     ! No. of sends for bcasts 
                                                  ! np_loop = numpro for MPI
      integer, allocatable, SAVE:: yfirst(:)      ! First latitude
      integer, allocatable, SAVE:: ylast(:)       ! Last latitude
      integer, allocatable, SAVE:: zfirst(:)      ! First level
      integer, allocatable, SAVE:: zlast(:)       ! Last level
      integer, SAVE:: ncall_r, ncall_s

!------------------------------------------------------------------------------
!  Variables to control global array locations and window synchronization
!------------------------------------------------------------------------------
      integer win_count                           ! Counts No. of windows in use
      integer lastwin                             ! ID of last synch'd window
      integer pkgs_per_pro                        ! No. of MPI packages per PE
      integer igosouth, igonorth                  ! Index of send direction
      integer ifromsouth, ifromnorth              ! Index of recv direction

!------------------------------------------------------------------------------
!  Local type declaration for mp_windows
!------------------------------------------------------------------------------
      type window
         integer :: id             ! Window id
         integer :: wsize          ! Size of global window (point based)
         integer :: ncall_s        ! Count send calls on window
         integer :: ncall_r        ! Count recv calls on window

#if defined(MPI2)
#if defined(LINUX)
#define MPI_OFFSET_KIND 8
#define MPI_ADDRESS_KIND 8
#endif
         integer(kind=MPI_ADDRESS_KIND) :: offset_s
         integer(kind=MPI_ADDRESS_KIND) :: offset_r
#else
         integer :: offset_s      ! Starting position in GA send
         integer :: offset_r      ! Starting position in GA recv
#endif
         integer :: dest          ! For use with send calls
         integer :: src           ! For use with recv calls
         integer :: size_r        ! Size of incoming message
     end type window

!------------------------------------------------------------------------------
! Beginning Global Array variable declaration:
!------------------------------------------------------------------------------

      type (window) :: r8_win
      type (window) :: r4_win
      type (window) :: i4_win

      integer ierror, irr

!-----------------------------
!  SHMEM variable declarations
!-----------------------------

#if defined(USE_SHMEM)
      integer ga_ptr
      integer, SAVE:: commglobal   ! Global Communicator

#define ga_t1_r ga_r8_r
#define ga_t1_s ga_r8_s
#define t1_win  r8_win

      real(r8), SAVE, TARGET:: ga_r8a_r(MAX( PLON*PLAT*(PLEV+1), idimsize*nbuf*max_call ))
      real(r8), SAVE, TARGET:: ga_r8a_s(MAX( PLON*PLAT*(PLEV+1), idimsize*nbuf*max_call ))

      real(r8), SAVE, TARGET:: ga_r8b_r(MAX( PLON*PLAT*(PLEV+1), idimsize*nbuf*max_call ))
      real(r8), SAVE, TARGET:: ga_r8b_s(MAX( PLON*PLAT*(PLEV+1), idimsize*nbuf*max_call ))

      real(r8), DIMENSION(:), SAVE, POINTER:: ga_r8_r
      real(r8), DIMENSION(:), SAVE, POINTER:: ga_r8_s

      real(r4), SAVE, TARGET:: ga_r4a_r(PLON*PLAT*PLEV)
      real(r4), SAVE, TARGET:: ga_r4a_s(PLON*PLAT*PLEV)

      real(r4), SAVE, TARGET:: ga_r4b_r(PLON*PLAT*PLEV)
      real(r4), SAVE, TARGET:: ga_r4b_s(PLON*PLAT*PLEV)

      real(r4), DIMENSION(:), SAVE, POINTER:: ga_r4_r
      real(r4), DIMENSION(:), SAVE, POINTER:: ga_r4_s

      integer(i4), SAVE, TARGET:: ga_i4a_r(PLON*PLAT*PLEV)
      integer(i4), SAVE, TARGET:: ga_i4a_s(PLON*PLAT*PLEV)

      integer(i4), SAVE:: ga_i4c_r(PLON*PLAT*PLEV)
      integer(i4), SAVE:: ga_i4c_s(PLON*PLAT*PLEV)

      integer(i4), SAVE, TARGET:: ga_i4b_r(PLON*PLAT*PLEV)
      integer(i4), SAVE, TARGET:: ga_i4b_s(PLON*PLAT*PLEV)

      integer(i4), DIMENSION(:), SAVE, POINTER:: ga_i4_r
      integer(i4), DIMENSION(:), SAVE, POINTER:: ga_i4_s

#else
# ifdef MPI2
!-----------------------------
!  MPI-2 variable declarations
!-----------------------------
      real(r8),    SAVE:: ga_t1_r(idimsize*nbuf*max_call)
      real(r8),    SAVE:: ga_t1_s(idimsize*nbuf*max_call)
      real(r8),    SAVE:: ga_r8_r(PLON*PLAT*PLEV)
      real(r8),    SAVE:: ga_r8_s(PLON*PLAT*PLEV)

      real(r4),    SAVE:: ga_r4_r(PLON*PLAT*PLEV)
      real(r4),    SAVE:: ga_r4_s(PLON*PLAT*PLEV)

      integer(i4), SAVE:: ga_i4_r(PLON*PLAT*PLEV)
      integer(i4), SAVE:: ga_i4_s(PLON*PLAT*PLEV)
      type (window) :: t1_win

#if defined(LINUX) || defined(OSF1)
#if defined(LINUX)
#define MPI_ADDRESS_KIND 8
#endif
      integer, parameter:: MPI_MODE_NOCHECK    = 0 
      integer, parameter:: MPI_MODE_NOSTORE    = 0    
      integer, parameter:: MPI_MODE_NOPUT      = 0 
      integer, parameter:: MPI_MODE_NOPRECEDE  = 0 
      integer, parameter:: MPI_MODE_NOSUCCEED  = 0 
#endif
      integer(kind=MPI_ADDRESS_KIND) bsize
      integer, SAVE:: commglobal   ! Global Communicator
      integer, SAVE:: Status(MPI_STATUS_SIZE)
#else
      real(r8), allocatable, SAVE:: ga_t1_r(:)
      real(r8), allocatable, SAVE:: ga_t1_s(:)
      real(r8), allocatable, SAVE:: ga_r8_r(:)
      real(r8), allocatable, SAVE:: ga_r8_s(:)

      real(r4), allocatable, SAVE:: ga_r4_r(:)
      real(r4), allocatable, SAVE:: ga_r4_s(:)

      integer(i4), allocatable, SAVE:: ga_i4_r(:)
      integer(i4), allocatable, SAVE:: ga_i4_s(:)
!-----------------------------
!  MPI-1 variable declarations
!-----------------------------
      type (window) :: t1_win
      integer, SAVE:: nsend                   ! Number of messages out-going
      integer, SAVE:: nrecv                   ! Number of messages in-coming
      integer, SAVE:: nread                   ! Number of messages read
      integer, allocatable, SAVE:: sqest(:)   ! Request handler for sends
      integer, allocatable, SAVE:: rqest(:)   ! Request handler for recvs
      integer, SAVE:: bsize             
      integer, SAVE:: commglobal              ! Global Communicator
      integer, SAVE:: Status(MPI_STATUS_SIZE)
      integer, allocatable, SAVE:: Stats(:)
#endif
#endif
!Balaji
      logical :: no_mod_comm=.TRUE. !TRUE on PEs that should NOT call mod_comm
                                    !set FALSE in mp_init by all PEs
                                    !that call mp_init

! !DESCRIPTION:
!
!    This module contains all SPMD parallelism decomposition and communication
!         routines.  MPI2 and MPI1 have been implemented. 
!
!    1) All implementations use real*8, real*4, and integer*4 windows
!       which are used with global arrays as follows:
!             r8_win -> ga_r8 - for use with real*8 types
!             r4_win -> ga_r4 - for use with real*4 types
!             i4_win -> ga_i4 - for use with integer*4 types
!
!       note: MPI routines need 2 buffers per GA, ga_<type>_s & ga_<type>_r
!             ga_<type>_r is used for the MPI2 windows
!
!    2) All global arrays are 1-dimensional, they are accessed
!       as needed inside the Ga_Put/Ga_Get routines with offset vars.
!       (Ga_Put/Ga_Get routines are all 4d with openmp on the 3rd (k) dim.)
!
! !REVISION HISTORY:
!    2001.09.01   Lin
!    2002.04.16   Putman  Modified for Global Array code
!    2002.04.16   Putman  Added ProTeX documentation
!    2002.05.28   Putman  Added use of precision module
!    2003.10.30   Lin     MLP removeid; dynamic memory allocation with MPI1
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_init --- Initialize SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_init( im, jm, km, comm )
!
! !INPUT PARAMETERS:
      integer, intent(in) :: im
      integer, intent(in) :: jm
      integer, intent(in) :: km
      integer, optional :: comm
! !DESCRIPTION:
!
!     Initialize SPMD parallel communication
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman	Modified for Global Array code
!    2002.04.09   Putman	Added ProTeX documentation
!    2002.08.06   Sawyer        Added optional communicator input argument                                               
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mysize

#if defined(USE_SHMEM)
        integer num_pes, my_pe
#endif


!Balaji: no_mod_comm is declared TRUE by default
!        but all PEs calling mp_init will set it to FALSE
! (i.e mp_init should be called only by those PEs that run mod_comm)
        no_mod_comm = .FALSE.
#if defined(MPI2) || defined(USE_SHMEM)
      if ( im*jm*km > PLON*PLAT*PLEV .or. im*(km+1) > PLON*(PLEV+1) ) then
         write(*,*) 'Warning in mp_init; buffer size decleared within mod_comm'
         write(*,*) 'is too small: blow up is to be expected'
      endif
#endif
        win_count = 0
 
#if defined(USE_SHMEM)
#if defined(IRIX64)
        CALL START_PES(0)
#else
        call SHMEM_INIT()
#endif

#if defined(USE_VT)
        call MPI_INIT(ierror)
#endif
        numpro = num_pes()
        gid = my_pe()
	ga_ptr = 1
	ga_r8_r => ga_r8a_r
	ga_r8_s => ga_r8a_s
	ga_r4_r => ga_r4a_r
	ga_r4_s => ga_r4a_s
	ga_i4_r => ga_i4a_r
	ga_i4_s => ga_i4a_s
	reduce_sync = SHMEM_SYNC_VALUE
	collect_sync = SHMEM_SYNC_VALUE
	bcast_sync = SHMEM_SYNC_VALUE
        call omp_start
        if( gid == 0 ) write(*,*) 'SHMEM is used within mod_comm'
#else
! MPI1 or MPI2

        if ( present(comm) ) then
          call mpi_start( comm )
        else
          call mpi_start( MPI_COMM_WORLD )
        endif

        call omp_start
#ifdef MPI2
        if( gid == 0 ) write(*,*) 'MPI2 is used within mod_comm'
#else
        if( gid == 0 ) write(*,*) 'MPI-1 is used within mod_comm'

        allocate( sqest(MAX(nbuf,numpro)*max_call) )
        allocate( rqest(MAX(nbuf,numpro)*max_call) )
        allocate( Stats(MAX(nbuf,numpro)*max_call*MPI_STATUS_SIZE) )

        plon  = im;  plat = jm;  plev = km
        idimsize = plon*nghost*(plev+1) 

        allocate ( ga_t1_r(idimsize*nbuf*max_call) )
        allocate ( ga_t1_s(idimsize*nbuf*max_call) )
        allocate ( ga_r8_r(PLON*PLAT*PLEV) )
        allocate ( ga_r8_s(PLON*PLAT*PLEV) )
        allocate ( ga_r4_r(PLON*PLAT*PLEV) )
        allocate ( ga_r4_s(PLON*PLAT*PLEV) )
        allocate ( ga_i4_r(PLON*PLAT*PLEV) )
        allocate ( ga_i4_s(PLON*PLAT*PLEV) )
#endif
        mysize = idimsize*nbuf*max_call
        call win_init_r8(t1_win, ga_t1_r, mysize)
#endif
        mysize = PLON*PLAT*PLEV
        call win_init_r8(r8_win, ga_r8_r, mysize)
        call win_init_r4(r4_win, ga_r4_r, mysize)
        call win_init_i4(i4_win, ga_i4_r, mysize)

        igosouth   = 0
        igonorth   = 1
        ifromsouth = 1
        ifromnorth = 0

        ncall_s = 0
        ncall_r = 0
        lastwin = r8_win%id

        np_loop = numpro

        allocate( yfirst( numpro ) )
        allocate( ylast ( numpro ) )
        allocate( zfirst( numpro ) )
        allocate( zlast ( numpro ) )
!EOC
      end subroutine mp_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_exit --- End SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_exit ( solo )

      logical, intent(in) :: solo          ! If mod_comm is the sole communication
                                           ! module or called within another framework
! !DESCRIPTION:
!
!     End SPMD parallel communication
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman	Modified for Global Array code
!    2002.04.09   Putman	Added ProTeX documentation
!    2003.11.18   Lin           Added input solo enabling mod_comm to be a sub-module
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if( no_mod_comm )return

        deallocate( yfirst )
        deallocate( ylast )
        deallocate( zfirst )
        deallocate( zlast )

#if !defined(USE_SHMEM)

#ifdef MPI2
        call MPI_WIN_FREE( t1_win%id, ierror )
        call MPI_WIN_FREE( r8_win%id, ierror )
        call MPI_WIN_FREE( r4_win%id, ierror )
        call MPI_WIN_FREE( i4_win%id, ierror )
#else
        deallocate( sqest )
        deallocate( rqest )
        deallocate( Stats )

        deallocate( ga_t1_r )
        deallocate( ga_t1_s )
        deallocate( ga_r8_r )
        deallocate( ga_r8_s )
        deallocate( ga_r4_r )
        deallocate( ga_r4_s )
        deallocate( ga_i4_r )
        deallocate( ga_i4_s )
#endif
        if(solo) call MPI_FINALIZE (ierror)

#elif defined(USE_SHMEM)
        call mp_barrier()
#if defined(USE_VT)
    	if(solo) call MPI_FINALIZE (ierror)
#endif
#endif
!EOC
      end subroutine mp_exit
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: omp_start --- Start openMP parallelism
!
! !INTERFACE:
      subroutine omp_start
! !DESCRIPTION:
!
!     Start openMP parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer ios, n, nowpro, nowcpu
        integer omp_get_num_threads
        character*80 evalue
        integer mp_suggested_numthreads

      if( no_mod_comm )return
#if defined(SET_CPUS)

!wfc moved up        character*80 evalue
        call getenv('AGCM_N_THREADS_PER_PROCESS',evalue)
        if (gid == 0) then
          read(evalue,*,iostat=ios) numcpu
          if ( ios .ne. 0 ) then
              print *, 'ERROR: cannot read AGCM_N_THREADS_PER_PROCESS', &
                       trim(evalue)
               call exit(1)
          end if
        endif
        call mp_bcst_int(numcpu)

#if defined(_OPENMP)

#if defined (IRIX64)
       call mp_set_numthreads(numcpu)  !keep it for a while, :)
#else
       call omp_set_num_threads(numcpu)
#endif

#else
       print *, 'ERROR: must turn on OPENMP to set numcpu in mod_comm'
       call exit(1)
#endif

#else 

#if defined(_OPENMP)

#if defined (IRIX64)
!wfc moved up        integer mp_suggested_numthreads
        numcpu = mp_suggested_numthreads(0)
#else
!$omp parallel
        numcpu = omp_get_num_threads()
!$omp end parallel
#endif

#else
        numcpu = 1
#endif

#endif

#if defined(MPI2) || defined(USE_SHMEM)
#if defined(MT_OFF)
        pkgs_per_pro = 1
#else
        pkgs_per_pro = numcpu
#endif
#endif

#if defined(PIN_CPUS)
       if (gid == 0 ) write(6,*) 'SGI specific CPU pinning activated!'
!$omp parallel do private(n,nowcpu)
        nowpro = gid
        do n=1,numcpu
          nowcpu = n + (nowpro) * numcpu-1
          call mp_assign_to_cpu(nowcpu)
        enddo
#endif
!EOC
      end subroutine omp_start
!------------------------------------------------------------------------------

#if !defined(USE_SHMEM)
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mpi_start --- Start MPI parallelism
!
! !INTERFACE:
      subroutine mpi_start( comm )
! !INPUT PARAMETERS:
      integer :: comm   ! Global communicator (may be MPI_COMM_WORLD)
! !DESCRIPTION:
!
!     Start MPI parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!    02.08.06   Sawyer  Added communicator input arguments
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        logical flag
        integer npthreads

      if( no_mod_comm )return
#if defined(MPI2) && !defined(AIX) && (!defined LINUX) && (!defined OSF1)
        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, npthreads, ierror)
        endif
        call MPI_QUERY_THREAD(npthreads, ierror)
#if !defined(MT_OFF)
        if (npthreads /= MPI_THREAD_MULTIPLE) then
!-- begin SJL ----
!         write(*,*) gid, 'did not provide MPI_THREAD_MULTIPLE. ', &
!               'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
!               'for multi-threading MPI2'
!
!         call MPI_FINALIZE(ierror)
!         call exit(1)
! -- End SJL mode ----
!         write(*,*) 'Calling MPI_INIT_THREAD ...'
          call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, npthreads, ierror)
          call MPI_QUERY_THREAD(npthreads, ierror)
!         write(*,*) 'npthreads=',npthreads, 'MPI_THREAD_MULTIPLE=',MPI_THREAD_MULTIPLE
        endif
#endif
#else
        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT( ierror )
        endif
#endif
        call MPI_COMM_RANK (comm, gid, ierror)
        call MPI_COMM_SIZE (comm, numpro, ierror)
        call MPI_COMM_DUP  (comm, commglobal, ierror)
!EOC
      end subroutine mpi_start
!------------------------------------------------------------------------------
#endif

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r8 --- Initialize real*8 communication window
!
! !INTERFACE:
      subroutine win_init_r8(win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: isize
        real(r8), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*8 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

      if( no_mod_comm )return
        win_count = win_count + 1

#if defined(MPI2)
        call MPI_INFO_CREATE(info, ierror)
        call MPI_INFO_SET(info, "no_locks", "true", ierror)
#if defined(AIX)
        info = MPI_INFO_NULL
#endif
        call MPI_TYPE_SIZE(mp_r8, mp_size, ierror)
        bsize = isize*mp_size
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, commglobal, &
                            win%id, ierror)
#if !defined(AIX)
        call MPI_INFO_FREE(info, ierror)
#endif
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#else
        win%id = win_count
#endif
        win%wsize = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r4 --- Initialize real*4 communication window
!
! !INTERFACE:
      subroutine win_init_r4(win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: isize
        real(r4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

      if( no_mod_comm )return
        win_count = win_count + 1

#if defined(MPI2)
        call MPI_INFO_CREATE(info, ierror)
        call MPI_INFO_SET(info, "no_locks", "true", ierror)
#if defined(AIX)
        info = MPI_INFO_NULL
#endif
        call MPI_TYPE_SIZE(mp_r4, mp_size, ierror)
        bsize = isize*mp_size
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, commglobal, &
                            win%id, ierror)
#if !defined(AIX)
        call MPI_INFO_FREE(info, ierror)
#endif
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#else
        win%id = win_count
#endif
        win%wsize = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_i4 --- Initialize integer*4 communication window
!
! !INTERFACE:
      subroutine win_init_i4(win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: isize
        integer(i4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize integer*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

      if( no_mod_comm )return
        win_count = win_count + 1

#if defined(MPI2)
        call MPI_INFO_CREATE(info, ierror)
        call MPI_INFO_SET(info, "no_locks", "true", ierror)
#if defined(AIX)
        info = MPI_INFO_NULL
#endif
        call MPI_TYPE_SIZE(mp_i4, mp_size, ierror)
        bsize = isize*mp_size
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, commglobal, &
                            win%id, ierror)
#if !defined(AIX)
        call MPI_INFO_FREE(info, ierror)
#endif
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#else
        win%id = win_count
#endif
        win%wsize = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_i4
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: y_decomp --- Decompose the Latitude & Level direction
!
! !INTERFACE:
      subroutine y_decomp(jm, km, jfirst, jlast, kfirst, klast, myid)
! !INPUT PARAMETERS:
      integer, intent(in):: jm     ! Dimensions
      integer, intent(in):: km     ! Levels
      integer, intent(in):: myid
! !OUTPUT PARAMETERS:
      integer, intent(out):: jfirst, jlast, kfirst, klast
! !DESCRIPTION:
!
!     Decompose the Latitude & Level direction for SPMD parallelism
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer p, p1, p2, lats, pleft
      integer, allocatable:: ydist(:)

      if( no_mod_comm )return
      if (myid == 0) then
          print *, "Two-level Parallelism: Domains x Threads=", numpro, numcpu
      endif

      allocate( ydist( numpro ) )

      lats = jm / numpro
      pleft = jm - lats * numpro


      do p=1,numpro
         ydist(p) = lats
      enddo

      if ( pleft .ne. 0 ) then
          p1 = (numpro+1) / 2 
          p2 = p1 + 1
        do while ( pleft .ne. 0 )
           if( p1 .eq. 1 ) p1 = numpro
               ydist(p1) = ydist(p1) + 1
               pleft = pleft - 1
               if ( pleft .ne. 0 ) then
                    ydist(p2) = ydist(p2) + 1
                    pleft = pleft - 1
               endif
               p2 = p2 + 1
               p1 = p1 - 1
        enddo
      endif

! Safety check:

      do p=1,numpro-1
      if( ydist(p) < 3 ) then
          write(*,*) 'Number of Proc is too large for jm=',jm
          call exit(1)
      endif
      enddo

      lats = 0
      do p = 1, numpro
         lats = lats + ydist(p)
      enddo

      if ( lats /= jm ) then
         print *, "Decomp: big trouble sum(ydist) = ", lats, "!=", jm
      endif
 
      jfirst = 1
      jlast  = ydist(1)
      yfirst(1) = jfirst
      ylast(1) = jlast
      kfirst = 1
      klast = km
      zfirst(1) = kfirst
      zlast(1) = klast

      do p = 1,numpro-1
         yfirst(p+1) = ylast(p) + 1
         ylast(p+1) = ylast(p) + ydist(p+1) 
         if( p == myid ) then
            jfirst = yfirst(p+1)
            jlast  = ylast (p+1)
         endif
         zfirst(p+1) = kfirst
         zlast(p+1) = klast
      enddo

      deallocate (ydist)
!EOC
      end subroutine y_decomp
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: set_decomp --- Set the Latitude & Level decomposition
!
! !INTERFACE:
      subroutine set_decomp(nprocs, jm, km, ydist, zdist)
! !INPUT PARAMETERS:
      integer, intent(in):: nprocs
      integer, intent(in):: jm, km
      integer, intent(in):: ydist(nprocs)
      integer, intent(in):: zdist(nprocs)   ! Currently not used
!
! !DESCRIPTION:
!
!     Set the Latitude & Level decomposition:
!             if it is defined external to mod_comm
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer lats, p

      if( no_mod_comm )return
! Safety check:
      lats = 0
      do p = 1, nprocs
         lats = lats + ydist(p)
      enddo

      if ( lats .ne. jm ) then
         print *, "Decomp: big trouble sum(ydist) = ", lats, "!=", jm
      endif

      yfirst(1) = 1
      ylast(1) = ydist(1)
      zfirst(1) = 1
      zlast(1) = km

      do p = 1,nprocs-1
         yfirst(p+1) = ylast(p) + 1
         ylast(p+1) = ylast(p) + ydist(p+1)
         zfirst(p+1) = 1
         zlast(p+1) = km
      enddo
!EOC
      end subroutine set_decomp
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send4d_ns --- Send 4d north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_send4d_ns(im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
      real(r8), intent(in):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Send 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(t1_win)

!--------------
! Send to south
!--------------

      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(MPI2)
        t1_win%src      = gid - 1
        t1_win%offset_r = ifromsouth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        t1_win%size_r   = im*ng_s*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest     = gid - 1
        t1_win%offset_s = igosouth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_Put4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst, jfirst+ng_n-1, kfirst, klast, 1, nq,   &
                         ga_t1_s, ga_t1_r )
      endif

!--------------
! Send to north
!--------------
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(MPI2)
        t1_win%src      = gid + 1
        t1_win%offset_r = ifromnorth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        t1_win%size_r   = im*ng_n*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest     = gid + 1
        t1_win%offset_s = igonorth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_Put4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast-ng_s+1, jlast, kfirst, klast, 1, nq,     &
                         ga_t1_s, ga_t1_r )
      endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_send4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv4d_ns --- Receive 4d north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_recv4d_ns(im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Receive 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Close(t1_win)


! Recv from south
      if ( jfirst > 1 ) then
        t1_win%src  = gid-1
        t1_win%offset_r = ifromsouth*idimsize + (t1_win%ncall_r-1)*idimsize*nbuf
        call Ga_Get4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        t1_win%src  = gid+1
        t1_win%offset_r = ifromnorth*idimsize + (t1_win%ncall_r-1)*idimsize*nbuf
        call Ga_Get4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast+1,     jlast+ng_n, kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_recv4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send2_ns --- Send 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_send2_ns(im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real(r8), intent(in):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(in):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 
!
! !DESCRIPTION:
!
!     Send 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(t1_win)

! Send to south
      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(MPI2)
        t1_win%src  = gid - 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromsouth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest = gid - 1
        t1_win%offset_s = igosouth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_Put4d_r8( q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(MPI2)
        t1_win%src  = gid + 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromnorth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest = gid + 1
        t1_win%offset_s = igonorth*idimsize + (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_Put4d_r8( q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jlast,     jlast,    kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jlast,     jlast,    kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_send2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv2_ns --- Receive 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_recv2_ns(im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(inout):: q2(im,jfirst-nd:jlast+nd,kfirst:klast)
!
! !DESCRIPTION:
!
!     Receive 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer j

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Close(t1_win)

! Recv from south
      if ( jfirst > 1 ) then
        j = jfirst - 1
        t1_win%src  = gid - 1
        t1_win%offset_r = ifromsouth*idimsize + (t1_win%ncall_r-1)*idimsize*nbuf
        call Ga_Get4d_r8( q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        j = jlast + 1
        t1_win%src  = gid + 1
        t1_win%offset_r = ifromnorth*idimsize + (t1_win%ncall_r-1)*idimsize*nbuf
        call Ga_Get4d_r8( q1, t1_win, im, jm, km, 2, & 
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_recv2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d --- Send ghost region
!
! !INTERFACE:
      subroutine mp_send3d(dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                  i1, i2, j1, j2, k1, k2, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(t1_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recv src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
      endif
#endif
! Send ghost region
      if ( dest >= 0 .and. dest < numpro ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_Put4d_r8( q, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
      endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_send3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d --- Recv ghost region
!
! !INTERFACE:
      subroutine mp_recv3d(src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                i1, i2, j1, j2, k1, k2, qout)
!
! !INPUT PARAMETERS:
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Close(t1_win)

! Recv from src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsize*nbuf
        call Ga_Get4d_r8( qout, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_recv3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d_2 --- Send 2 ghost regions
!
! !INTERFACE:
      subroutine mp_send3d_2(dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                        i1, i2, j1, j2, k1, k2, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q1(if:il, jf:jl, kf:kl)
      real(r8), intent(in):: q2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send two general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(t1_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recv src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + t1_win%size_r 
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
      endif
#endif
! Send ghost region
      if ( dest >= 0 .and. dest < numpro ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsize*nbuf
        call Ga_Put4d_r8( q1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Put4d_r8( q2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,  &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_s, ga_t1_r  )
      endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_send3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d_2 --- Recv 2 ghost regions
!
! !INTERFACE:
      subroutine mp_recv3d_2(src, im, jm, km, if, il, jf, jl, kf, kl, &
                                  i1, i2, j1, j2, k1, k2, qout1, qout2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout1(if:il, jf:jl, kf:kl)
      real(r8), intent(inout):: qout2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv two general 3d real*8 ghost regions
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return

#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Close(t1_win)

! Recv from src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsize*nbuf
        call Ga_Get4d_r8( qout1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Get4d_r8( qout2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,   &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_recv3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_reduce_sum --- Sum real*8 scalar over all processes
!
! !INTERFACE:
      subroutine mp_reduce_sum (sum)
!
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: sum
!
! !DESCRIPTION:
!
!     Sum real*8 scalar over all processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer  n
      real(r8) sumin
      real(r8) sumg(numpro)

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined(USE_SHMEM)
      call Win_Open(r8_win)
      ga_r8_s(1) = sum
      call SHMEM_REAL8_SUM_TO_ALL (ga_r8_r, ga_r8_s, 1, 0, 0, numpro, &
                             ga_r8_s(2), reduce_sync)
      call Win_Close(r8_win)
      sum = ga_r8_r(1)
      call Win_Finalize(r8_win)
#else
      sumin = sum
      call MPI_ALLREDUCE( sumin, sum, 1, mp_r8, MPI_SUM, &
                          commglobal, ierror )
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_reduce_sum
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_reduce_max --- Reduce maximum real*8 scalar over all processes
!
! !INTERFACE:
      subroutine mp_reduce_max(km, cymax)
!-----
! !INPUT PARAMETERS:
      integer, intent(in):: km
!
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: cymax(km)
!
! !DESCRIPTION:
!
!     Reduce maximum real*8 scalar over all processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer k, n
      real(r8) maxin(km)

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined(USE_SHMEM)
      call Win_Open(r8_win)
      do k=1,km
         ga_r8_s(k) = cymax(k)
      enddo
      call SHMEM_REAL8_MAX_TO_ALL (ga_r8_r, ga_r8_s, km, 0, 0, numpro, &
                                   ga_r8_s(km+1), reduce_sync)
      call Win_Close(r8_win)
      do k=1,km
         cymax(k) = ga_r8_r(k)
      enddo
      call Win_Finalize(r8_win)
#else
!$omp parallel do private(k, irr)
      do k=1,km
        maxin(k) = cymax(k)
      enddo
      call MPI_ALLREDUCE( maxin, cymax, km, mp_r8, MPI_MAX, &
                          commglobal, irr )
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_reduce_max
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_minmax --- Reduce min/max real*8 scalar over all processes
!
! !INTERFACE:
      subroutine mp_minmax(qmin, qmax)
!
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qmin, qmax
!
! !DESCRIPTION:
!
!     Reduce minimum/maximum real*8 scalar over all processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real(r8) minin, maxin
      integer n

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined(USE_SHMEM)
      call Win_Open(r8_win)
      ga_r8_s(1) = qmax
      ga_r8_s(2) = -1.*qmin
      call SHMEM_REAL8_MAX_TO_ALL (ga_r8_r, ga_r8_s, 2, 0, 0, numpro, &
                             ga_r8_s(3), reduce_sync)
      call Win_Close(r8_win)
      qmax = ga_r8_r(1)
      qmin = -1.*ga_r8_r(2)
      call Win_Finalize(r8_win)
#else
      maxin = qmax
      minin = qmin
      call MPI_ALLREDUCE(maxin, qmax, 1, mp_r8, MPI_MAX, &
                         commglobal, ierror)
      call MPI_ALLREDUCE(minin, qmin, 1, mp_r8, MPI_MIN, &
                         commglobal, ierror)
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_minmax
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sum1d --- Sum real*8 over all latitudes
!
! !INTERFACE:
      subroutine mp_sum1d(jm, jfirst, jlast, qin, sum0)
!
! !INPUT PARAMETERS:
      integer, intent(in):: jm
      integer, intent(in):: jfirst, jlast
      real(r8), intent(in):: qin(jfirst:jlast)
! !OUTPUT PARAMETERS:
      real(r8), intent(out)::  sum0
!
! !DESCRIPTION:
!
!     Sum real*8 over all latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer j, n
      real(r8) qout(jm)

      if( no_mod_comm )return
#if defined(QUICK_DIRTY_CODE)
     sum0 = 0.
     do j=jfirst,jlast
        sum0 = sum0 + qin(j)
     enddo
     call mp_reduce_sum(sum0)
#else
      call mp_gather3d(qin, qout, 1, jm, 1, jfirst, jlast, &
                       1, 1, 0, 0, 0)
! Compute the sum if "Master"
      if ( gid == 0 ) then
         sum0 = 0.
         do j=1,jm
            sum0 = sum0 + qout(j)
         enddo
      endif
      call mp_bcst_real(sum0)
#endif
!EOC
      end subroutine mp_sum1d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_add1d --- Sum real*8 array over all processes
!
! !INTERFACE:
      subroutine mp_add1d(jdim, qin)
!
! !INPUT PARAMETERS:
      integer, intent(in):: jdim
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qin(jdim)
!
! !DESCRIPTION:
!
!     Sum real*8 array over all processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real(r8) qout(jdim,numpro)
      integer j, n

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(r8_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      if (gid==0) then
        do n=1,numpro
          r8_win%src = n-1
          r8_win%size_r = jdim 
          r8_win%offset_r = r8_win%src*jdim
          call Ga_RecvInit_r8(r8_win, ga_r8_r)
        enddo
      endif
#endif

      r8_win%dest = 0
      r8_win%offset_s = gid*jdim
      call Ga_Put4d_r8(qin, r8_win, 1, 1, jdim, numpro, &
                            1, 1, 1, 1, 1, jdim, gid+1, gid+1, &
                            1, 1, 1, 1, 1, jdim, gid+1, gid+1, &
                            ga_r8_s, ga_r8_r  )

      call Win_Close(r8_win)

      if (gid==0) then
        do n=1,numpro
          r8_win%src = n-1
          r8_win%offset_r = r8_win%src*jdim
          call Ga_Get4d_r8(qout(1,n), r8_win, 1, 1, jdim, numpro, &
                         1, 1, 1, 1, 1, jdim, r8_win%src+1, r8_win%src+1, &
                         1, 1, 1, 1, 1, jdim, r8_win%src+1, r8_win%src+1, &
                         ga_r8_r  )
        enddo
        do j=1,jdim
          qin(j) = 0.
          do n=1,numpro
            qin(j) = qin(j) + qout(j,n)
          enddo
        enddo
      endif

      call Win_Finalize(r8_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_add1d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_bcst_real --- Broadcast real*8 scalar
!
! !INTERFACE:
      subroutine mp_bcst_real(val)
!
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: val
!
! !DESCRIPTION:
!
!     Broadcast real*8 scalar from Process=id to All other Processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
      real(r8) val_arr(1)

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

     val_arr(1) = val
     call Ga_Broadcast_r8(val_arr, 1)
     val = val_arr(1)

!EOC
      end subroutine mp_bcst_real
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_bcst_int --- Broadcast integer*4 scalar
!
! !INTERFACE:
      subroutine mp_bcst_int(intv)
!
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: intv
!
! !DESCRIPTION:
!
!     Broadcast integer*4 scalar from Process=id to All other Processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
       integer n
       integer(i4) intv_arr(1)

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif
   
      intv_arr(1) = intv
      call Ga_Broadcast_i4(intv_arr, 1)
      intv = intv_arr(1)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_bcst_int
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_bcst_r2d --- Broadcast real*8 2d array
!
! !INTERFACE:
      subroutine mp_bcst_r2d(im, jm, jfirst, jlast, qin, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm
      integer, intent(in):: id        ! source ID
      integer, intent(in):: jfirst, jlast
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qin(im,jm)
!
! !DESCRIPTION:
!
!     Broadcast real*8 2d array from Process=id to All other Processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Ga_Broadcast_r8 ( qin, im*jm ) 

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_bcst_r2d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_gath_r2d --- Gather real*8 2d array
!
! !INTERFACE:
      subroutine mp_gath_r2d(im, jm, jfirst, jlast, qin, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm
      integer, intent(in):: id        ! source ID
      integer, intent(in):: jfirst, jlast
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qin(im,jm)
!
! !DESCRIPTION:
!
!     Gather real*8 2d array from All Processes to Process=id
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, j1, j2
 
      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(r8_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      if (gid == id) then
         do n=1,numpro
           r8_win%src  = n-1
           r8_win%size_r = im*(ylast(n)-yfirst(n)+1)
           r8_win%offset_r = (yfirst(n)-1)*im
           call Ga_RecvInit_r8(r8_win, ga_r8_r)
         enddo
      endif
#endif

      r8_win%dest = id
      r8_win%offset_s = (jfirst-1)*im
      call Ga_Put4d_r8 ( qin, r8_win, 1, im, jm, 1, &
                         1, 1, 1, im, 1,      jm,    1, 1, &
                         1, 1, 1, im, jfirst, jlast, 1, 1, &
                         ga_r8_s, ga_r8_r  )

      call Win_Close(r8_win)

      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           r8_win%src  = n-1
           r8_win%offset_r = (j1-1)*im
           call Ga_Get4d_r8 ( qin, r8_win, 1, im, jm, 1, &
                              1, 1, 1, im, 1,  jm, 1, 1, &
                              1, 1, 1, im, j1, j2, 1, 1, &
                              ga_r8_r  )
         enddo
      endif

      call Win_Finalize(r8_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_gath_r2d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_gath_3d_int --- Gather integer*4 3d array
!
! !INTERFACE:
      subroutine mp_gath_3d_int(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: idim, jdim, kdim
      integer, intent(in):: id
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout)::    q(idim,jdim,kdim)
!
! !DESCRIPTION:
!
!     Gather integer*4 3d array from All Processes to Process=id,
!            with non-standard (ie, /= PLON, PLAT, PLEV) array dimensions.
!
! !REVISION HISTORY:
!    2002.09.09   Putman        
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, ndims
      parameter (ndims=6)
      integer(i4) all_dims(ndims, numpro)
      integer ir1, ir2, jr1, jr2, kr1, kr2

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      all_dims(1,gid+1) = i1
      all_dims(2,gid+1) = i2
      all_dims(3,gid+1) = j1
      all_dims(4,gid+1) = j2
      all_dims(5,gid+1) = k1
      all_dims(6,gid+1) = k2

      call Ga_AllToAll_i4(all_dims, numpro*ndims, ndims, gid*ndims+1)

      call Win_Open(i4_win)
#if !defined(USE_SHMEM) && !defined(MPI2)
      if (gid==id) then
        i4_win%offset_r = 0
        do n=1,numpro
          i4_win%src = n-1
          ir1 = all_dims(1,n)
          ir2 = all_dims(2,n)
          jr1 = all_dims(3,n)
          jr2 = all_dims(4,n)
          kr1 = all_dims(5,n)
          kr2 = all_dims(6,n)
          if (n>1) then
            i4_win%offset_r = i4_win%offset_r + &
                             ((all_dims(2,n-1)-all_dims(1,n-1)+1) * &
                              (all_dims(4,n-1)-all_dims(3,n-1)+1) * &
                              (all_dims(6,n-1)-all_dims(5,n-1)+1))
          endif
          i4_win%size_r = (ir2-ir1+1)*(jr2-jr1+1)*(kr2-kr1+1)
          call Ga_RecvInit_i4(i4_win, ga_i4_r)
        enddo
      endif
#endif
      i4_win%dest = id
      i4_win%offset_s = 0
      do n=1,gid
        i4_win%offset_s = i4_win%offset_s + &
                         ((all_dims(2,n)-all_dims(1,n)+1) * &
                          (all_dims(4,n)-all_dims(3,n)+1) * &
                          (all_dims(6,n)-all_dims(5,n)+1))
      enddo
      call Ga_Put4d_i4(q, i4_win, idim, jdim, kdim, 1, &
                          1,  idim, 1,  jdim, 1,  kdim, 1, 1, &
                          i1, i2,   j1, j2,   k1, k2,   1, 1, &
                          ga_i4_s, ga_i4_r  )
      call Win_Close(i4_win)
      if (gid==id) then
        i4_win%offset_r = 0
        do n=1,numpro
          i4_win%src = n-1
          ir1 = all_dims(1,n)
          ir2 = all_dims(2,n)
          jr1 = all_dims(3,n)
          jr2 = all_dims(4,n)
          kr1 = all_dims(5,n)
          kr2 = all_dims(6,n)
          if (n>1) then
            i4_win%offset_r = i4_win%offset_r + &
                             ((all_dims(2,n-1)-all_dims(1,n-1)+1) * &
                              (all_dims(4,n-1)-all_dims(3,n-1)+1) * &
                              (all_dims(6,n-1)-all_dims(5,n-1)+1))
          endif
          call Ga_Get4d_i4(q, i4_win, idim, jdim, kdim, 1, &
                              1,   idim, 1,   jdim, 1,   kdim, 1, 1, &
                              ir1, ir2,  jr1, jr2,  kr1, kr2,  1, 1, &
                              ga_i4_r  )
        enddo
      endif
      call Win_Finalize(i4_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_gath_3d_int
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_gath_3d --- Gather real*8 3d array
!
! !INTERFACE:
      subroutine mp_gath_3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: idim, jdim, kdim
      integer, intent(in):: id
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout)::    q(idim,jdim,kdim)
!
! !DESCRIPTION:
!
!     Gather real*8 3d array from All Processes to Process=id, 
!            with non-standard (ie, /= PLON, PLAT, PLEV) array dimensions.
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, ndims
      parameter (ndims=6)
      integer(i4) all_dims(ndims, numpro)
      integer ir1, ir2, jr1, jr2, kr1, kr2

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      all_dims(1,gid+1) = i1
      all_dims(2,gid+1) = i2
      all_dims(3,gid+1) = j1
      all_dims(4,gid+1) = j2
      all_dims(5,gid+1) = k1
      all_dims(6,gid+1) = k2

      call Ga_AllToAll_i4(all_dims, numpro*ndims, ndims, gid*ndims+1)

      call Win_Open(r8_win)
#if !defined(USE_SHMEM) && !defined(MPI2)
      if (gid==id) then
        r8_win%offset_r = 0
        do n=1,numpro
          r8_win%src = n-1
          ir1 = all_dims(1,n)
          ir2 = all_dims(2,n)
          jr1 = all_dims(3,n)
          jr2 = all_dims(4,n)
          kr1 = all_dims(5,n)
          kr2 = all_dims(6,n)
          if (n>1) then
            r8_win%offset_r = r8_win%offset_r + & 
                             ((all_dims(2,n-1)-all_dims(1,n-1)+1) * &
                              (all_dims(4,n-1)-all_dims(3,n-1)+1) * &
                              (all_dims(6,n-1)-all_dims(5,n-1)+1))
          endif
          r8_win%size_r = (ir2-ir1+1)*(jr2-jr1+1)*(kr2-kr1+1)
          call Ga_RecvInit_r8(r8_win, ga_r8_r)
        enddo
      endif
#endif
      r8_win%dest = id
      r8_win%offset_s = 0
      do n=1,gid
        r8_win%offset_s = r8_win%offset_s + &
                         ((all_dims(2,n)-all_dims(1,n)+1) * &
                          (all_dims(4,n)-all_dims(3,n)+1) * &
                          (all_dims(6,n)-all_dims(5,n)+1))
      enddo

      call Ga_Put4d_r8(q, r8_win, idim, jdim, kdim, 1, &
                          1,  idim, 1,  jdim, 1,  kdim, 1, 1, &
                          i1, i2,   j1, j2,   k1, k2,   1, 1, &
                          ga_r8_s, ga_r8_r  )
      call Win_Close(r8_win)
      if (gid==id) then
        r8_win%offset_r = 0
        do n=1,numpro
          r8_win%src = n-1
          ir1 = all_dims(1,n)
          ir2 = all_dims(2,n)
          jr1 = all_dims(3,n)
          jr2 = all_dims(4,n)
          kr1 = all_dims(5,n)
          kr2 = all_dims(6,n)
          if (n>1) then
            r8_win%offset_r = r8_win%offset_r + &
                             ((all_dims(2,n-1)-all_dims(1,n-1)+1) * &
                              (all_dims(4,n-1)-all_dims(3,n-1)+1) * &
                              (all_dims(6,n-1)-all_dims(5,n-1)+1))
          endif
          call Ga_Get4d_r8(q, r8_win, idim, jdim, kdim, 1, &
                              1,   idim, 1,   jdim, 1,   kdim, 1, 1, &
                              ir1, ir2,  jr1, jr2,  kr1, kr2,  1, 1, &
                              ga_r8_r  )
        enddo
      endif
      call Win_Finalize(r8_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_gath_3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_scat3d_int --- Scatter integer*4 3d array
!
! !INTERFACE:
      subroutine mp_scat3d_int(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: idim, jdim, kdim
      integer, intent(in):: id
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: q(idim,jdim,kdim)
!
! !DESCRIPTION:
!
!     Scatter integer*4 3d array from Process=id to All other Processes, 
!            with non-standard (ie, /= PLON, PLAT, PLEV) array dimensions.
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i, j, k, n

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(i4_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      i4_win%src  = id
      i4_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      i4_win%offset_r = 0
      call Ga_RecvInit_i4(i4_win, ga_i4_r)
#endif

      if (gid == id) then
        do n=1,numpro
          i4_win%dest = n-1
          i4_win%offset_s = 0
          call Ga_Put4d_i4 ( q, i4_win, idim, jdim, kdim, 1, &
                               1,  idim, 1,  jdim, 1,  kdim, 1, 1, &
                               i1, i2,   j1, j2,   k1, k2,   1, 1, &
                               ga_i4_s, ga_i4_r  )
        enddo
      endif

      call Win_Close(i4_win)

      i4_win%src  = id
      i4_win%offset_r = 0
      call Ga_Get4d_i4 ( q, i4_win, idim, jdim, kdim, 1, &
                           1,  idim, 1,  jdim, 1,  kdim, 1, 1, &
                           i1, i2,   j1, j2,   k1, k2,   1, 1, &
                           ga_i4_r  )

      call Win_Finalize(i4_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_scat3d_int
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_scat3d --- Scatter real*8 3d array
!
! !INTERFACE:
      subroutine mp_scat3d(q, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: idim, jdim, kdim
      integer, intent(in):: id
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout)::    q(idim,jdim,kdim)
!
! !DESCRIPTION:
!
!     Scatter real*8 3d array from Process=id to All other Processes, 
!            with non-standard (ie, /= PLON, PLAT, PLEV) array dimensions.
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i, j, k, n

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(r8_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      r8_win%src  = id
      r8_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      r8_win%offset_r = 0
      call Ga_RecvInit_r8(r8_win, ga_r8_r)
#endif

      if (gid == id) then
        do n=1,numpro
          r8_win%dest = n-1
          r8_win%offset_s = 0
          call Ga_Put4d_r8 ( q, r8_win, idim, jdim, kdim, 1, &
                                1,  idim, 1,  jdim, 1,  kdim, 1, 1, &
                                i1, i2,   j1, j2,   k1, k2,   1, 1, &
                                ga_r8_s, ga_r8_r  )
        enddo
      endif

      call Win_Close(r8_win)

      r8_win%src  = id
      r8_win%offset_r = 0
      call Ga_Get4d_r8 ( q, r8_win, idim, jdim, kdim, 1, &
                            1,  idim, 1,  jdim, 1,  kdim, 1, 1, &
                            i1, i2,   j1, j2,   k1, k2,   1, 1, &
                            ga_r8_r  )

      call Win_Finalize(r8_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_scat3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_gath4_r2d --- Gather real*4 2d array
!
! !INTERFACE:
      subroutine mp_gath4_r2d(im, jm, jfirst, jlast, qin, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm
      integer, intent(in):: id        ! source ID
      integer, intent(in):: jfirst, jlast
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: qin(im,jm)
!
! !DESCRIPTION:
!
!     Gather real*4 2d array from Process=id to All other Processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, j1, j2

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(r4_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      if (gid == id) then
         do n=1,numpro
           r4_win%src  = n-1
           r4_win%size_r = im*(ylast(n)-yfirst(n)+1) 
           r4_win%offset_r = (yfirst(n)-1)*im
           call Ga_RecvInit_r4(r4_win, ga_r4_r)
         enddo
      endif
#if defined(OSF1)
      call mp_barrier
#endif
#endif

      r4_win%dest = id
      r4_win%offset_s = (jfirst-1)*im
      call Ga_Put4d_r4 ( qin, r4_win, 1, im, jm, 1, &
                         1, 1, 1, im, 1,      jm,    1, 1, &
                         1, 1, 1, im, jfirst, jlast, 1, 1, &
                         ga_r4_s, ga_r4_r  )

      call Win_Close(r4_win)

      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           r4_win%src  = n-1
           r4_win%offset_r = (j1-1)*im
           call Ga_Get4d_r4 ( qin,  r4_win, 1, im, jm, 1, &
                              1, 1, 1, im, 1,  jm, 1, 1, &
                              1, 1, 1, im, j1, j2, 1, 1, &
                              ga_r4_r  )
         enddo
      endif

      call Win_Finalize(r4_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_gath4_r2d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_barrier --- Synchronize all SPMD processes
!
! !INTERFACE:
      subroutine mp_barrier
!
! !DESCRIPTION:
!
!     Synchronize all SPMD processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if( no_mod_comm )return
#if defined (USE_SHMEM)
        call SHMEM_BARRIER_ALL()
#else
        call MPI_BARRIER(commglobal, ierror)
#endif
!EOC
      end subroutine mp_barrier
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_scatter3d --- Scatter real*8 4d array
!
! !INTERFACE:
      subroutine mp_scatter3d(qin, qout, im, jm, km, jfirst, jlast, &
                              kfirst, klast, ng_s, ng_n, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: ng_s, ng_n
      integer, intent(in):: id           ! Source (usually 0)
      real(r8), intent(in):: qin(im,jm,km)
! !OUTPUT PARAMETERS:
      real(r8), intent(out):: qout(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
!
! !DESCRIPTION:
!
!     Scatter real*8 4d array from Process=id to All other Processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2001.11.15   Putman	Modified for Global Arrays code
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n
      integer j1, j2, k1, k2
 
      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(r8_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      r8_win%src = id
      r8_win%size_r = im*(jlast-jfirst+1)*(klast-kfirst+1)
      r8_win%offset_r = ((jfirst-1)*km+(kfirst-1))*im
      call Ga_RecvInit_r8(r8_win, ga_r8_r)
#endif

      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           if (km == PLEV) then
             k1 = zfirst(n)
             k2 = zlast(n)
           else
             k1 = 1
             k2 = km
           endif
           r8_win%dest = n-1
           r8_win%offset_s = ((j1-1)*km+(k1-1))*im
           call Ga_Put4d_r8 ( qin, r8_win, im, jm, km, 1, &
                              1, im, 1,  jm, 1,  km, 1,1, &
                              1, im, j1, j2, k1, k2, 1,1, &
                              ga_r8_s, ga_r8_r  )
         enddo
      endif

      call Win_Close(r8_win)

      r8_win%offset_r = ((jfirst-1)*km+(kfirst-1))*im
      r8_win%src = id
      call Ga_Get4d_r8 ( qout, r8_win, im, jm, km, 1, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, 1, &
                         1, im, jfirst,      jlast,      kfirst, klast, 1, 1, &
                         ga_r8_r  )

      call Win_Finalize(r8_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_scatter3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_gather3d --- Gather real*8 4d array
!
! !INTERFACE:
      subroutine mp_gather3d(qin, qout, im, jm, km, jfirst, jlast, &
                             kfirst, klast, ng_s, ng_n, id)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: ng_s, ng_n
      integer, intent(in):: id         ! process ID to gather data to
      real(r8), intent(in):: qin(im,jfirst-ng_s:jlast+ng_n,kfirst:klast)
! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(im,jm,km)
!
! !DESCRIPTION:
!
!     Gather real*8 4d array from Process=id to All other Processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2001.11.15   Putman	Modified for Global Arrays code
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i, j, k, iq
      integer j1, j2, k1, k2
      integer n

      if( no_mod_comm )return
#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

      call Win_Open(r8_win)

#if !defined(USE_SHMEM) && !defined(MPI2)
! Init Recvs
      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           if (km == PLEV) then
             k1 = zfirst(n)
             k2 = zlast(n)
           else
             k1 = 1
             k2 = km
           endif
           r8_win%src  = n-1
           r8_win%size_r = im*(j2-j1+1)*(k2-k1+1)
           r8_win%offset_r = ((j1-1)*km+(k1-1))*im
           call Ga_RecvInit_r8(r8_win, ga_r8_r)
         enddo
      endif
#endif

      r8_win%dest = id
      r8_win%offset_s = ((jfirst-1)*km+(kfirst-1))*im
      call Ga_Put4d_r8 ( qin, r8_win, im, jm, km,  1, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1,  1, &
                         1, im, jfirst,      jlast,      kfirst, klast, 1,  1, &
                         ga_r8_s, ga_r8_r  )

      call Win_Close(r8_win)

      if (gid == id) then
         do n=1,numpro
           j1 = yfirst(n)
           j2 = ylast(n)
           if (km == PLEV) then
             k1 = zfirst(n)
             k2 = zlast(n)
           else
             k1 = 1
             k2 = km
           endif
           r8_win%src  = n-1
           r8_win%offset_r = ((j1-1)*km+(k1-1))*im
           call Ga_Get4d_r8 ( qout, r8_win, im, jm, km,  1, &
                              1, im,  1, jm,  1, km, 1,  1, &
                              1, im, j1, j2, k1, k2, 1,  1, &
                              ga_r8_r  )
         enddo
      endif

      call Win_Finalize(r8_win)

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_gather3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_wrt3d --- Write real*8 3d array to file
!
! !INTERFACE:
      subroutine mp_wrt3d(iout, nrec, r_zero, qin, im, jm, km, &
                          jfirst, jlast, kfirst, klast, id )
!
! !INPUT PARAMETERS:
      integer, intent(in):: iout
      integer, intent(in):: nrec
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast, kfirst, klast
      integer, intent(in):: id         ! process ID
      real(r8), intent(in):: r_zero
      real(r8), intent(in):: qin(im,jfirst:jlast,kfirst:klast)
!
! !DESCRIPTION:
!
!     Write real*8 3d array to 4-byte unit number iout, record=nrec:
!           If < r_zero write 0.
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2001.11.15   Putman	Modified for Global Arrays code
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i, j, k, iq
      real(r4) a2(im,jm)
      integer n
      integer qsize
      integer send_tag, recv_tag
      real(r8) qout(im,jm,km)

      real(r8) qmax


      if( no_mod_comm )return
      call mp_gather3d(qin, qout, im, jm, km, jfirst, jlast, &
                       kfirst, klast, 0, 0, 0)
      if (gid == id) then
        do k=1,km
          do j=1,jm
            do i=1,im
              if( abs(qout(i,j,k)) < r_zero ) then
                 qout(i,j,k) = 0.
              endif
              a2(i,j) = qout(i,j,k)
            enddo
          enddo
        write(iout,rec=nrec+k) a2
        enddo
      endif
!EOC
      end subroutine mp_wrt3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Open --- Open a communication window
!
! !INTERFACE:
      subroutine Win_Open(win)
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Begin a communication epoch, by opening a comm window.
!     Update number of send calls on the window (win%ncall_s).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
      win%ncall_s = win%ncall_s + 1
      ncall_s = ncall_s + 1

#if defined(USE_SHMEM)
      if (ncall_s == 1) then
         if (ga_ptr == 1) then 
            ga_ptr = 2
            ga_r8_r => ga_r8b_r
            ga_r8_s => ga_r8b_s
            ga_r8_r => ga_r8b_r
            ga_r8_s => ga_r8b_s
            ga_r4_r => ga_r4b_r
            ga_r4_s => ga_r4b_s
            ga_i4_r => ga_i4b_r
            ga_i4_s => ga_i4b_s
         else
            ga_ptr = 1
            ga_r8_r => ga_r8a_r
            ga_r8_s => ga_r8a_s
            ga_r8_r => ga_r8a_r
            ga_r8_s => ga_r8a_s
            ga_r4_r => ga_r4a_r
            ga_r4_s => ga_r4a_s
            ga_i4_r => ga_i4a_r
            ga_i4_s => ga_i4a_s
         endif
      endif
#endif

!EOC
      end subroutine Win_Open
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Close --- Close a communication window
!
! !INTERFACE:
      subroutine Win_Close(win)
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     End a communication epoch, by closing a comm window.
!     Update number of receive calls on the window (win%ncall_r).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
      win%ncall_r = win%ncall_r + 1
      ncall_r = ncall_r + 1
#if defined(USE_SHMEM)
      if (ncall_r == 1) then
          call mp_barrier()
      endif
#endif
#if defined(MPI2)
      if (win%ncall_r == 1) then
          call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                             win%id, ierror)
      endif
#endif

!EOC
      end subroutine Win_Close
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Finalize --- Reset a communication window after a comm epoch.
!
! !INTERFACE:
      subroutine Win_Finalize(win)
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Complete a communication epoch and reset a comm window.
!     Update global lastwin with win%id.
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY:
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if( no_mod_comm )return
      if (ncall_s == ncall_r) then
#if !defined(USE_SHMEM) && !defined(MPI2)
        call MPI_WAITALL(nsend, sqest, Stats, ierror)
        nsend = 0
        nrecv = 0
        nread = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif

      if (win%ncall_s == win%ncall_r) then
#if defined(MPI2)
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#endif
        lastwin = win%id
        win%ncall_s = 0
        win%ncall_r = 0
      endif

!EOC
      end subroutine Win_Finalize
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r8 --- Write to real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r8 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2

! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga_s(win%wsize)
      real(r8), intent(inout):: ga_r(win%wsize)
!
! !DESCRIPTION:
!
!     Write to real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
      integer send_tag, qsize

#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
#if defined(OSF1)
      integer(kind=8) p, tmpsize, mysize, mydisp
#else
      integer p, tmpsize, mysize, mydisp
#endif
#endif

      if( no_mod_comm )return
      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length 
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2 and MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                 ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)

#if defined(MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))

!$omp parallel do private(p,mysize,mydisp, irr)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize

#if defined(MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_r8, &
                     win%dest, mydisp, mysize, mp_r8, &
                     win%id, irr)
#else
           call SHMEM_PUT64(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
#endif
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r8, win%dest, &
                     send_tag, commglobal, sqest(nsend), irr)
#endif

!EOC
      end subroutine Ga_Put4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r8 --- Initiate real*8 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r8( win, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga(win%wsize)
!
! !DESCRIPTION:
!
!     Initiate real*8 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

      if( no_mod_comm )return
#if !defined(USE_SHMEM) && !defined(MPI2)
      if (win%wsize >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        nrecv    = nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r8, win%src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
      endif
#endif

!EOC
      end subroutine Ga_RecvInit_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r8 --- Read from real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r8 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: ga(win%wsize)
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      if( no_mod_comm )return
#if !defined(USE_SHMEM) && !defined(MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r4 --- Write to real*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r4 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: ga_s(win%wsize)
      real(r4), intent(inout):: ga_r(win%wsize)
!
! !DESCRIPTION:
!
!     Write to real*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
      integer send_tag, qsize
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
#if defined(OSF1)
      integer(kind=8) p, tmpsize, mysize, mydisp
#else
      integer p, tmpsize, mysize, mydisp
#endif
#endif

      if( no_mod_comm )return
      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2 and MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))

!$omp parallel do private(p,mysize,mydisp,irr)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
#if defined(MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_r4, &
                     win%dest, mydisp, mysize, mp_r4, &
                     win%id, irr)
#else
           call SHMEM_PUT32(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
#endif
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r4, win%dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r4 --- Initiate real*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r4( win, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: ga(win%wsize)
!
! !DESCRIPTION:
!
!     Initiate real*4 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

      if( no_mod_comm )return
#if !defined(USE_SHMEM) && !defined(MPI2)
      recv_tag = win%src
      qsize    = win%size_r
      nrecv    = nrecv + 1
      call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r4, win%src, &
                     recv_tag, commglobal, rqest(nrecv), ierror)
#endif
!EOC
      end subroutine Ga_RecvInit_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r4 --- Read from real*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r4 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r4), intent(in)  :: ga(win%wsize)
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
      if( no_mod_comm )return
#if !defined(USE_SHMEM) && !defined(MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_i4 --- Write to integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_i4 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga_s(win%wsize)
      integer(i4), intent(inout):: ga_r(win%wsize)
!
! !DESCRIPTION:
!
!     Write to integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
      integer send_tag, qsize
#if defined(MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
#if defined(OSF1)
      integer(kind=8) p, tmpsize, mysize, mydisp
#else
      integer p, tmpsize, mysize, mydisp
#endif
#endif

      if( no_mod_comm )return
      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2 and MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,mydisp, irr)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
#if defined(MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_i4, &
                     win%dest, mydisp, mysize, mp_i4, &
                     win%id, irr)
#else
           call SHMEM_PUT32(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
#endif
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_i4, win%dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_i4 --- Initiate integer*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_i4( win, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga(win%wsize)
!
! !DESCRIPTION:
!
!     Initiate integer*4 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

      if( no_mod_comm )return
#if !defined(USE_SHMEM) && !defined(MPI2)
      recv_tag = win%src
      qsize    = win%size_r
      nrecv    = nrecv + 1
      call MPI_IRECV(ga(win%offset_r+1), qsize, mp_i4, win%src, &
                     recv_tag, commglobal, rqest(nrecv), ierror)
#endif
!EOC
      end subroutine Ga_RecvInit_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_i4 --- Read from integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_i4 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: ga(win%wsize)
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
      if( no_mod_comm )return
#if !defined(USE_SHMEM) && !defined(MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r8 --- Broadcast an real*8 1d global array
!
! !INTERFACE:
        subroutine Ga_Broadcast_r8 ( q, isize )
! !INPUT PARAMETERS:
        integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
        real(r8), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!	Broadcast an real*8 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
        integer i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) p, tmpsize, mysize, mydisp
#else
        integer p, tmpsize, mysize, mydisp
#endif
#endif

      if( no_mod_comm )return
#if defined(USE_SHMEM)
        mysize = isize
	call mp_barrier
        do i=1,mysize
           ga_r8_s(i) = q(i)
        enddo
	call SHMEM_BROADCAST64(ga_r8_r, ga_r8_s, mysize, 0, 0, 0, numpro, bcast_sync)
	if (gid /= 0) then
	do i=1,isize
           q(i) = ga_r8_r(i)
	enddo
	endif
#endif
#if !defined(USE_SHMEM)
      call MPI_BCAST(q, isize, mp_r8, 0, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r4 --- Broadcast an real*4 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_r4 ( q, isize )
! !INPUT PARAMETERS:
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!	Broadcast an real*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
        integer i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) mysize
#else
        integer mysize
#endif
#endif

      if( no_mod_comm )return
#if defined(USE_SHMEM)
	mysize = isize
	call mp_barrier
        do i=1,mysize
           ga_r4_s(i) = q(i)
        enddo
	call SHMEM_BROADCAST32(ga_r4_r, ga_r4_s, mysize, 0, 0, 0, numpro, bcast_sync)
	if (gid /= 0) then
	do i=1,isize
           q(i) = ga_r4_r(i)
	enddo
	endif
#endif
#if !defined(USE_SHMEM)
      call MPI_BCAST(q, isize, mp_r4, 0, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_i4 --- Broadcast an integer*4 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_i4 ( q, isize )
! !INPUT PARAMETERS:
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!	Broadcast an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
        integer i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) mysize
#else
        integer mysize
#endif
#endif

      if( no_mod_comm )return
#if defined(USE_SHMEM)
	mysize = isize
	call mp_barrier
        do i=1,mysize
           ga_i4_s(i) = q(i)
        enddo
	call SHMEM_BROADCAST32(ga_i4_r, ga_i4_s, mysize, 0, 0, 0, numpro, bcast_sync)
	if (gid /= 0) then
	do i=1,isize
           q(i) = ga_i4_r(i)
	enddo
	endif
#endif
#if !defined(USE_SHMEM)
      call MPI_BCAST(q, isize, mp_i4, 0, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r8 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
        subroutine Ga_AllToAll_r8 ( q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
        integer, intent(in)  :: Gsize   ! Global size of array
        integer, intent(in)  :: Lsize   ! size of Local portion
        integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
        real(r8), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!	All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
        integer n, i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) my_Lsize
#else
        integer my_Lsize
#endif
#endif

      if( no_mod_comm )return
#if defined(USE_SHMEM)
	my_Lsize = Lsize
	call mp_barrier
        do i=1,Lsize
           ga_r8_s(i) = q(i+istart-1)
        enddo
	call SHMEM_COLLECT8(ga_r8_r, ga_r8_s, my_Lsize, 0, 0, numpro, collect_sync)
	do i=1,Gsize
	   q(i) = ga_r8_r(i)
	enddo
#endif
#if !defined(USE_SHMEM)
        call MPI_ALLGATHER(q(istart), Lsize, mp_r8, q, Lsize, mp_r8, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_AllToAll_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r4 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_r4 ( q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
        integer, intent(in)  :: Gsize   ! Global size of array
        integer, intent(in)  :: Lsize   ! size of Local portion
        integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
        real(r4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!	All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
         integer n, i
#if defined(USE_SHMEM)
#if defined(OSF1)
         integer(kind=8) my_Lsize
#else
         integer my_Lsize
#endif
#endif

      if( no_mod_comm )return
#if defined(USE_SHMEM)
	my_Lsize = Lsize
	call mp_barrier
        do i=1,Lsize
           ga_r4_s(i) = q(i+istart-1)
        enddo
	call SHMEM_COLLECT4(ga_r4_r, ga_r4_s, my_Lsize, 0, 0, numpro, collect_sync)
	do i=1,Gsize
	   q(i) = ga_r4_r(i)
	enddo
#endif
#if !defined(USE_SHMEM)
        call MPI_ALLGATHER(q(istart), Lsize, mp_r4, q, Lsize, mp_r4, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_AllToAll_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_i4 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
        subroutine Ga_AllToAll_i4( q, Gsize, Lsize, istart)
! !INPUT PARAMETERS:
        integer, intent(in)  :: Gsize   ! Global size of array
        integer, intent(in)  :: Lsize   ! size of Local portion
        integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
        integer(i4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!	All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
        integer n, i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) my_Lsize
#else
        integer my_Lsize
#endif
#endif

      if( no_mod_comm )return
#if defined(USE_SHMEM)
	my_Lsize = Lsize
	call mp_barrier
        do i=1,Lsize
           ga_i4_s(i) = q(i+istart-1)
        enddo
	call SHMEM_COLLECT4(ga_i4_r, ga_i4_s, my_Lsize, 0, 0, numpro, collect_sync)
	do i=1,Gsize
	   q(i) = ga_i4_r(i)
	enddo
#endif
#if !defined(USE_SHMEM)
        call MPI_ALLGATHER(q(istart), Lsize, mp_i4, q, Lsize, mp_i4, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_AllToAll_i4
!------------------------------------------------------------------------------


#if defined(USE_MPI_IO)
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_file_open --- Open a binary file for writing with MPI-IO
!
! !INTERFACE:
      subroutine mp_file_open ( unit, filename )
! !INPUT PARAMETERS:
      character(*), intent(in) :: filename
! !OUTPUT PARAMETERS:
      integer, intent(inout) :: unit
!
! !DESCRIPTION:
!
!     Open a binary file for writing with MPI-IO
!
! !REVISION HISTORY:
!    2002.06.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      if( no_mod_comm )return
      call MPI_FILE_OPEN( commglobal, trim(filename),  &
                          MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                          MPI_INFO_NULL, unit, ierror )

!EOC
      end subroutine mp_file_open
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_file_write_at --- Write real*4 data to a binary file with MPI-IO
!
! !INTERFACE:
      subroutine mp_file_write_at ( unit, nrec, g_size, l_size, g_start, qout )
! !INPUT PARAMETERS:
      integer, intent(in) :: unit    ! File unit number
      integer, intent(in) :: nrec    ! Current record number
      integer, intent(in) :: g_size  ! Size of global array being writen
      integer, intent(in) :: l_size  ! Local size of global array this proc is writing
      integer, intent(in) :: g_start ! location in global array this proc writes to
      real(r4), intent(in):: qout(l_size)  ! Local data to be written
!
! !DESCRIPTION:
!
!     Write real*4 data to a binary file with MPI-IO
!
! !REVISION HISTORY:
!    2002.06.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer(MPI_OFFSET_KIND) offset 
      integer io_stat(MPI_STATUS_SIZE)

      if( no_mod_comm )return
      offset = nrec-1
      call MPI_FILE_WRITE_AT( unit, (offset*g_size + g_start)*4, &
                              qout, l_size, MPI_REAL, io_stat, ierror )

!EOC
      end subroutine mp_file_write_at
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_file_close --- CLose a binary file for writing with MPI-IO
!
! !INTERFACE:
      subroutine mp_file_close ( unit )
! !INPUT PARAMETERS:
      integer, intent(in) :: unit
!
! !DESCRIPTION:
!
!     Close a binary file for writing with MPI-IO
!
! !REVISION HISTORY:
!    2002.06.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      if( no_mod_comm )return
      call MPI_FILE_CLOSE( unit, ierror )

!EOC
      end subroutine mp_file_close
!------------------------------------------------------------------------------
#endif

#endif

end module mod_comm


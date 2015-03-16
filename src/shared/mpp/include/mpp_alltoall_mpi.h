subroutine MPP_ALLTOALL_(sbuf, rbuf, pelist)

    MPP_TYPE_, dimension(:), intent(in) :: sbuf
    MPP_TYPE_, dimension(:), intent(inout) :: rbuf

    integer, intent(in), optional :: pelist(0:)
    integer :: n

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    n = get_peset(pelist)
    if (peset(n)%count .eq. 1) return

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    if (verbose) call mpp_error(NOTE, 'MPP_ALLTOALL_: using MPI_Alltoall...')

    call MPI_Alltoall(sbuf, 1, MPI_TYPE_, rbuf, 1, MPI_TYPE_, &
                      peset(n)%id, error)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALL_

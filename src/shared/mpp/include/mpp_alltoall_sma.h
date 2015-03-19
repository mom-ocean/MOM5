subroutine MPP_ALLTOALL_(sbuf, rbuf, pelist)

    MPP_TYPE_, dimension(:), intent(in) :: sbuf
    MPP_TYPE_, dimension(:), intent(inout) :: rbuf

    integer, intent(in), optional :: pelist(0:)
    integer :: n

    call mpp_error(FATAL, 'MPP_ALLTOALL: No SHMEM implementation.')

    return
end subroutine MPP_ALLTOALL_

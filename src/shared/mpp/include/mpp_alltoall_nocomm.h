subroutine MPP_ALLTOALL_(sbuf, rbuf, pelist)

    MPP_TYPE_, dimension(:), intent(in) :: sbuf
    MPP_TYPE_, dimension(:), intent(inout) :: rbuf

    integer, intent(in), optional :: pelist(0:)
    integer :: n

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    return
end subroutine MPP_ALLTOALL_

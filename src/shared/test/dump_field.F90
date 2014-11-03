
module dump_field

! Dump a field. Subsequent calls with the same field name and proc will a new
! time point to the same file. See the test program below for use. 

use netcdf
implicit none

private 
public dump_field_2d, dump_field_close

type field_info_type
    integer :: count
    integer :: ncid
    integer :: varid
    integer :: min_varid
    integer :: max_varid
    integer :: mean_varid
    character(len=256) :: field_name
end type field_info_type

! Array containing list of field names and associated data, e.g. ncid.
integer, parameter :: MAX_FIELDS = 256
type (field_info_type), dimension(MAX_FIELDS) :: field_info
integer :: field_num = 1

contains  

subroutine setup_2d(field_name, proc_num, nx, ny, do_full_dump)
    
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: nx, ny, proc_num
    logical, intent(in) :: do_full_dump

    integer :: ncid, varid, min_varid, max_varid, mean_varid
    integer :: x_dimid, y_dimid, t_dimid
    integer :: status
    character(len=6) :: proc_str

    if (proc_num > 999999) then
        stop 'Error: dump_field::setup_2d(), proc_num too high.'
    endif

    write(proc_str, '(I6.6)') proc_num

    ! Open a file, set up a meta-data needed to save the field. 
    status = nf90_create(field_name//'.'//proc_str//'.nc', NF90_CLOBBER, ncid)
    call check(status, ' create file for '//field_name)
    call check(nf90_def_dim(ncid, 't', NF90_UNLIMITED, t_dimid), &
                            ' define t dim for '//field_name)
    if (do_full_dump) then 
        call check(nf90_def_dim(ncid, 'x', nx, x_dimid), &
                                'define x dim for '//field_name)
        call check(nf90_def_dim(ncid, 'y', ny, y_dimid), &
                                'define y dim for'//field_name)
        call check(nf90_def_var(ncid, field_name, NF90_DOUBLE, &
                                (/ x_dimid, y_dimid, t_dimid /), varid), &
                                ' define var for '//field_name)
        field_info(field_num)%varid = varid
    endif
    call check(nf90_def_var(ncid, field_name//'_min', NF90_DOUBLE, &
                            (/ t_dimid /), min_varid), &
                            ' define min var for '//field_name)
    call check(nf90_def_var(ncid, field_name//'_max', NF90_DOUBLE, &
                            (/ t_dimid /), max_varid), &
                            ' define max var for '//field_name)
    call check(nf90_def_var(ncid, field_name//'_mean', NF90_DOUBLE, &
                            (/ t_dimid /), mean_varid), &
                            ' define mean var for '//field_name)
    call check(nf90_enddef(ncid), ' enddef for '//field_name)

    field_info(field_num)%field_name = field_name
    field_info(field_num)%ncid = ncid 
    field_info(field_num)%min_varid = min_varid
    field_info(field_num)%max_varid = max_varid
    field_info(field_num)%mean_varid = mean_varid
    field_info(field_num)%count = 1

    field_num = field_num + 1
    if (field_num > MAX_FIELDS) then
        stop 'Error: dump_field::setup_2d(), too many fields.'
    endif

end subroutine

subroutine dump_field_2d(field_name, proc_num, field_data, do_full_dump)
    ! Search through field names and find appropriate file handle.
    ! Save to file. 

    character(len=*), intent(in) :: field_name
    integer, intent(in) :: proc_num
    real, dimension(:,:), intent(in) :: field_data
    logical, intent(in), optional :: do_full_dump

    real :: mean, divisor
    integer :: start(3), data_size(3), idx
    integer :: status
    logical :: found, dump

    found = .false.
    dump = .false.

    if (present(do_full_dump)) then 
        dump = do_full_dump
    end if

    call get_index(field_name, idx, found)
    if (.not. found) then
        call setup_2d(field_name, proc_num, size(field_data, 1), &
                      size(field_data, 2), dump)
    end if

    call get_index(field_name, idx, found)
    if (.not. found) then
        stop 'dump_field_mod::field_write'
    end if

    data_size = (/ size(field_data, 1), size(field_data, 2), 1 /)
    start = (/ 1, 1, field_info(idx)%count /)

    ! Dump data
    if (dump) then 
        status = nf90_put_var(field_info(idx)%ncid, field_info(idx)%varid, &
                              field_data, start=start, count=data_size)
        call check(status, ' put var for '//field_name)
    end if

    ! Write out some stats. 
    status = nf90_put_var(field_info(idx)%ncid, field_info(idx)%max_varid, &
                          (/ maxval(field_data) /), &
                          start=(/ field_info(idx)%count /), count=(/ 1 /))
    call check(status, ' put max var for '//field_name)
           
    status = nf90_put_var(field_info(idx)%ncid, field_info(idx)%min_varid, &
                          (/ minval(field_data) /), &
                          start=(/ field_info(idx)%count /), count=(/ 1 /))
    call check(status, ' put min var for '//field_name)
    divisor = size(field_data, 1) * size(field_data, 2)
    if (divisor /= 0) then 
        mean = sum(field_data) / divisor
    else
        mean = 0.0
    endif
    status = nf90_put_var(field_info(idx)%ncid, field_info(idx)%mean_varid, &
                         (/ mean /), &
                         start=(/ field_info(idx)%count /), count=(/ 1 /))
    call check(status, ' put mean var for '//field_name)
               
    call check(nf90_sync(field_info(idx)%ncid), ' sync file for '//field_name)

    field_info(idx)%count = field_info(idx)%count + 1 

end subroutine

subroutine dump_field_close(field_name)

    character(len=*), intent(in) :: field_name

    integer :: idx
    logical :: found

    call get_index(field_name, idx, found)
    if (.not. found) then
        stop 'dump_field_mod::field_close()'
    end if

    call check(nf90_close(field_info(idx)%ncid))
    
end subroutine

subroutine get_index(field_name, idx, found)

    character(len=*), intent(in) :: field_name
    integer, intent(out) :: idx
    logical, intent(out) :: found

    found = .false.
    do idx=1, field_num
        if (field_name == field_info(idx)%field_name) then
            found = .true.
            return 
        end if
    end do
   
end subroutine get_index

subroutine check(status, msg)
    integer, intent ( in) :: status
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: error_msg

    error_msg = 'dump_field_mod::check() '//nf90_strerror(status)
    if (present(msg)) then 
        error_msg = trim(error_msg)//' at: '//msg
    end if

    if(status /= nf90_noerr) then 
        print *, error_msg
        stop 'dump_field_mod::check()'
    end if
end subroutine check  

end module

!program test_dump_field

!   use dump_field

!   real, dimension(3, 3) :: array

!   array = reshape((/ 12.0, 32.0, 1.23123, 55.0, 3322.0, 65.0, 0.123, 99.0, 10.0 /), shape(array))
!   call dump_field_2d('sst', 4, array, .true.)

!   array = reshape((/ 3099.0, 554343.0, 0.3221, 405.0, 23.0, 56.0, 123.0, 89.0, 0.000002 /), shape(array))
!   call dump_field_2d('sst', 4, array, .true.)

!   call dump_field_close('sst')

!end program


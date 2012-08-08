module check_mask_mod

  use mpp_mod,         only : stdout, mpp_error, WARNING, FATAL, mpp_pe, mpp_root_pe, get_unit
  use mpp_domains_mod, only : mpp_define_layout, mpp_compute_extent
  use fms_mod,         only : open_namelist_file, close_file, check_nml_error, file_exist, string

  implicit none 
  private

  !--- namelist interface
  integer :: max_pe = 128 
  integer :: halo = 1
  logical :: show_valid_only = .false.
  namelist /check_mask_nml/ max_pe, halo, show_valid_only
!--- namelist variables for OBC check----------------------------------------------

  integer, parameter :: max_obc = 4      ! maximum number of open boundaries (increase if want more)
  integer                               :: nobc= 0
  character(len=10), dimension(max_obc) :: direction='none' ! open directions; to be consistent with is, ie, js, je
  integer, dimension(max_obc)           :: is=-999, ie=-999, js=-999, je=-999     ! boundary position
  integer, dimension(max_obc)           :: nsmooth=2        ! number of points to be smoothed
  character(len=32), dimension(max_obc) :: name='none' 
  
  namelist /obc_nml/ nobc, direction, is, ie, js, je, nsmooth, name


  public :: checking_mask

contains

  subroutine checking_mask (wet_in, cyclic_x, cyclic_y, is_tripolar, have_obc)   
    real,     intent(in) :: wet_in(:,:)
    logical,  intent(in) :: cyclic_x, cyclic_y, is_tripolar
    logical,  intent(in) :: have_obc

    integer              :: unit, ierr, io, nx, ny
    integer, allocatable :: mask(:,:), mask_list(:,:), obc_error(:,:)
    integer, allocatable :: ibegin(:), iend(:), jbegin(:), jend(:)
    real,    allocatable :: wet(:,:)
    integer              :: isd, ied, jsd, jed
    integer              :: isc, iec, jsc, jec
    integer              :: isg, ieg, jsg, jeg
    integer              :: i, j, layout(2), no, np, nmask, n, ii, jj
    character(len=64)    :: formatstring
    logical              :: on_domain = .false.
    character(len=128)   :: tablefile
    integer              :: fileunit, io_stat

    if(file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, nml=check_mask_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'check_mask_nml')  ! also initializes nml error codes
       enddo
10     call close_file (unit)
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, nml=obc_nml, iostat=io, end=20)
          ierr = check_nml_error(io,'obc_nml')  ! also initializes nml error codes
       enddo
20     call close_file (unit)
    endif

    layout = 0
    nx = size(wet_in,1)
    ny = size(wet_in,2)

    allocate(wet(1-halo:nx+halo,1-halo:ny+halo) )
    wet = 0.0
    wet(1:nx,1:ny) = wet_in(1:nx,1:ny)

    ! fill the global halo.
    if(cyclic_x) then
       wet(1-halo:0, 1:ny) = wet(nx-halo+1:nx, 1:ny)
       wet(nx+1:nx+halo, 1:ny) = wet(1:halo, 1:ny)
    end if

    if(cyclic_y) then
       wet(1:nx, 1-halo:0) = wet(1:nx,ny-halo+1:ny)
       wet(1:nx, ny+1:ny+halo) = wet(1:nx, 1:halo)
    end if

    if(is_tripolar) wet(1-halo:nx+halo,ny+1:ny+halo) = wet(nx+halo:1-halo:-1,ny:ny-halo+1:-1)

    write(stdout(),'(a)') 'Checking for possible masking:'
    write(stdout(),*)     'Assuming ',halo,' halo rows. '
    write(stdout(),*)     'Total domain size: ', nx, ny
    if(show_valid_only .and. have_obc) then
      write(stdout(),*)     'Only layouts valid with OBC are shown.'
    else
      if (have_obc) write(stdout(),*)     'Also layouts invalid with OBC are shown.'
    endif

    isg = 1 
    ieg = nx 
    jsg = 1 
    jeg = ny 

    do np = 4, min(max_pe,nx*ny)
       call mpp_define_layout((/1,nx,1,ny/),np,layout)
       if (layout(1) > nx .or. layout(2) > ny) cycle
       allocate ( mask(layout(1), layout(2)) )
       if (have_obc) allocate ( obc_error(layout(1), layout(2)) )
       allocate(ibegin(layout(1)), iend(layout(1)) )
       allocate(jbegin(layout(2)), jend(layout(2)) )
       call mpp_compute_extent(1,nx,layout(1),ibegin,iend)
       call mpp_compute_extent(1,ny,layout(2),jbegin,jend)

       mask   = 0
       if (have_obc) obc_error = 0
       do j = 1, layout(2)
          jsc = jbegin(j) 
          jec = jend(j) 
          jsd = jbegin(j) - halo
          jed = jend(j) + halo
          do i = 1, layout(1)
             isc = ibegin(i) 
             iec = iend(i) 
             isd = ibegin(i) - halo
             ied = iend(i) + halo
             if (ANY(wet(isd:ied,jsd:jed) > 0.5) ) mask(i,j) = 1
	     if (have_obc) then
	       do no=1, nobc
! check, if the open boundary is at a domain     
                 on_domain=.false.
		 do jj = jsd, jed
                   do ii = isd, ied
                    if(jj.ge.js(no).and.jj.le.je(no).and.ii.ge.is(no).and.ii.le.ie(no)) on_domain = .true.
                   enddo
                 enddo
		 if (.not.on_domain) cycle
! check, if the open boundary is at a domain halo	     
                 select case( trim(direction(no)) )
                 case ('west') 
                   if (isg /= isc) then
                     if ((is(no)-isc) < 0)    obc_error(i,j)=no
                   endif
                   if ((ied-is(no)) <= halo)  obc_error(i,j)=no
                 case ('east') 
                   if (ieg /= iec) then
                     if ((iec-is(no)) < 0)    obc_error(i,j)=no
                   endif
                   if ((is(n)-isd) <= halo)   obc_error(i,j)=no
                 case ('south') 
                   if (jsg /= jsc) then
                     if ((js(no)-jsc) < 0)    obc_error(i,j)=no
                   endif
                   if ((jec-js(no)) <= halo)  obc_error(i,j)=no
                 case ('north')
		   if (jeg /= jec) then
                     if ((jec-js(no)) < 0)    obc_error(i,j)=no
                   endif
                   if ((js(no)-jsc) <= halo)  obc_error(i,j)=no
                 end select
	       enddo
	     endif
          enddo
       enddo

       nmask = count(mask == 0)
       if(nmask > 0) then
          allocate(mask_list(2, nmask))
          n = 0
          do j = 1, layout(2)
             do i = 1, layout(1)
                if(mask(i,j) == 0) then
                   n = n + 1
                   mask_list(1,n) = i
                   mask_list(2,n) = j
                end if
             end do
          end do

          if(nmask == np) call mpp_error(WARNING, "check_mask_mod: all points are land")
          if(show_valid_only .and. have_obc) then 
	    if(ANY(obc_error > 0)) then
              deallocate(ibegin, iend, jbegin, jend)
	      deallocate(mask, obc_error)
              deallocate(mask_list)
	      cycle
	    endif    
	  endif    

          ! write the n_mask, layout_mask and mask_list information to a table file
          if(mpp_pe() == mpp_root_pe()) then
             tablefile = "mask_table."//trim(string(nmask))//"."//trim(string(layout(1)))//"x"//trim(string(layout(2)))
             fileunit = get_unit()
             open(unit=fileunit, file=trim(tablefile), iostat=io_stat)
             if(io_stat .NE. 0) call mpp_error(FATAL,"check_mask: error in opening file"//trim(tablefile) )
             write(fileunit, *)trim(string(nmask))
             write(fileunit,*)trim(string(layout(1))), " , ", trim(string(layout(2)))
             do n = 1, nmask
                write(fileunit,*) trim(string(mask_list(1,n))), ",", trim(string(mask_list(2,n)))
             enddo
             close(unit=fileunit, iostat=io_stat)
             if(io_stat .NE. 0) call mpp_error(FATAL,"check_mask: error in closing file"//trim(tablefile) )
          endif

          write(stdout(),*) '_______________________________________________________________________'
          write(stdout(),*) 'Possible setting to mask out all-land points region, for use in coupler_nml'
          write(stdout(),'((a),i4)') 'Total number of domains                                  ', np 
          write(stdout(),'((a),i4)') 'Number of tasks (excluded all-land region) to be used is ', np - nmask
          write(stdout(),'((a),i4)') 'Number of regions to be masked out                       ', nmask
          write(stdout(),'((a),2i4)') 'The layout is                                           ', layout
!          write(stdout(),*) 'The list of masked region mask_list is ( specified x-layout and y-layout)'
!          write(stdout(),*) 'mask_list = ', mask_list
          write(stdout(),*) 'Masked and used tasks, 1: used, 0: masked'
	  write(formatstring,*) "(",layout(2),"(2x,",layout(1),"i1,/))"
          write(stdout(),formatstring) ((mask(i,j),i=1,layout(1)),j=layout(2),1,-1)
          if (have_obc .and. ANY(obc_error > 0)) then
	     write(stdout(),*) 'OBC at halos is not allowed yet. Do not use this configuration.'
	     write(stdout(),formatstring) ((obc_error(i,j),i=1,layout(1)),j=layout(2),1,-1)
          endif
	  write(*,*) ' domain decomposition'
          write (*,110) (iend(i)-ibegin(i)+1, i= 1, layout(1))
          write (*,120) (jend(i)-jbegin(i)+1, i= 1, layout(2))
110       format ('  X-AXIS = ',24i4,/,(11x,24i4))
120       format ('  Y-AXIS = ',24i4,/,(11x,24i4))
          write(stdout(),*) 'used, masked, layout ', np - nmask, nmask, layout 
          write(stdout(),*) ' '
          write(stdout(),*) 'To chose this mask layout please put the following lines in coupler_nml '
          write(stdout(),130) nmask
          write(stdout(),140) layout(1), layout(2)
          if (nmask > 1) then
	    write(formatstring,*) "('     mask_list = ',i0,',',i0,",nmask-1,"(',',i0,',',i0))"
          else
	    write(formatstring,*) "('     mask_list = ',i0,',',i0)"
	  endif
          write(stdout(),trim(formatstring)) ((mask_list(i,j),i=1,2),j=1,nmask)
130       format ('     n_mask = ',i4)
140       format ('     layout_mask = ',i4,',',i4)
          deallocate(mask_list)
       end if
       deallocate(ibegin, iend, jbegin, jend)
       deallocate (mask)
       if (have_obc) deallocate (obc_error)
    enddo

    return

  end subroutine checking_mask


end module check_mask_mod

! program to check mask for existing grid file.

#ifdef check_mask
program check_mask_driver
  use mpp_mod,        only : stdout
  use mpp_io_mod,     only : mpp_open,  MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
  use mpp_io_mod,     only : mpp_get_info, atttype, mpp_get_atts, mpp_close
  use mpp_io_mod,     only : mpp_get_att_name, mpp_get_att_char
  use fms_mod,        only : open_namelist_file, close_file, check_nml_error, file_exist
  use fms_mod,        only : fms_init, fms_end, field_size, read_data
  use check_mask_mod, only : checking_mask
implicit none

  integer :: dims(4), unit, ierr, io, ndim, nvar, natt, ntime, i, ni, nj
  logical :: cyclic_x, cyclic_y, is_tripolar
  type(atttype), allocatable :: global_atts(:)
  real,          allocatable :: wet(:,:)
  logical                    :: have_obc = .false.
  character(len=128) :: grid_file = "INPUT/grid_spec.nc"
  namelist /check_mask_driver_nml/ grid_file, have_obc


  call fms_init()

  if(file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr=1
     do while (ierr /= 0)
        read  (unit, nml=check_mask_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'check_mask_driver_nml')  ! also initializes nml error codes
     enddo
10   call close_file (unit)
  endif

  ! get boundary condition

  call mpp_open(unit,trim(grid_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
  call mpp_get_info(unit, ndim, nvar, natt, ntime)

  cyclic_x = .FALSE.
  cyclic_y = .FALSE.
  is_tripolar = .FALSE.
  allocate(global_atts(natt))
  call mpp_get_atts(unit,global_atts)
  do i=1,natt
     select case (trim(mpp_get_att_name(global_atts(i))))
     case ('x_boundary_type')
        if (trim(mpp_get_att_char(global_atts(i))) == 'cyclic') cyclic_x = .true.
     case ('y_boundary_type')
        if (trim(mpp_get_att_char(global_atts(i))) == 'fold_north_edge') then
           is_tripolar = .true.
        else if(trim(mpp_get_att_char(global_atts(i))) == 'cyclic')  then
           cyclic_y = .true.
        end if
     end select
  end do

  deallocate(global_atts)
  call mpp_close(unit)

  ! read land/sea mask information from grid file.
  call field_size(grid_file, 'wet', dims )
  ni = dims(1)
  nj = dims(2)
  allocate(wet(ni,nj))
  call read_data(grid_file, 'wet', wet  )
  call checking_mask(wet, cyclic_x, cyclic_y, is_tripolar, have_obc) ! The namelist for obc is read in checking_mask
  write(stdout(),*) " successfully run program check_mask"

  call fms_end()

end program check_mask_driver

#endif

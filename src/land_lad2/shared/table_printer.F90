module table_printer_mod
implicit none ; private

! TODO: add alignment to the list of arguments
! TODO: make possible to print line numbers

! ---- public interface
public :: table_printer_type
public :: init_with_headers  ! initializes the table with specified headers
public :: add_row
public :: print
public :: dealloc

! ---- module constants   
character(*), parameter :: DEFAULT_FORMAT      = '(g16.9)'
integer     , parameter :: DEFAULT_CELL_WIDTH  = 16
integer     , parameter :: DEFAULT_TABLE_WIDTH = HUGE(1)
integer     , parameter :: N_ROWS_INCREMENT    = 16

! ---- module types
type :: table_printer_type
   private
   integer   :: nr = 0
   integer   :: nc = 0
   character(32), pointer :: table(:,:) => null() ! table cells
end type

! ---- generic procedure interface
interface add_row
   module procedure add_row_real, add_row_integer
   module procedure add_row_logical, add_row_character
end interface add_row

contains 


! ==============================================================================
! deallocats all table data and resets the counters
subroutine dealloc(t)
   type(table_printer_type), intent(inout) :: t
   
   if (associated(t%table))  deallocate(t%table)
   t%table => null() ; t%nr = 0 ; t%nc = 0
end subroutine 


! ==============================================================================
! sets the column headers and resets the table data
subroutine init_with_headers(t,names)
   type(table_printer_type), intent(inout) :: t
   character(*), intent(in) :: names(:)
   
   call dealloc(t)
   t%nc = size(names)
   call realloc_table(t)
   t%table(1:,0) = names(:)
end subroutine 


! ==============================================================================
! increases a number of row storage. Private procedure.
subroutine realloc_table(t)
   type(table_printer_type), intent(inout) :: t

   character(32), pointer :: ptr(:,:)
   integer :: i,j,n

   ! reallocate data
   n = 0 ; if (associated(t%table)) n = ubound(t%table,2)
   if (t%nr >= n) then
      allocate(ptr(0:t%nc,0:t%nr+N_ROWS_INCREMENT))
      forall (i=lbound(ptr,1):ubound(ptr,1), j=lbound(ptr,2):ubound(ptr,2)) ptr(i,j) = ''
      if (associated(t%table)) then
         ptr(:,0:t%nr) = t%table(:,0:t%nr)
         deallocate(t%table)
      endif
      t%table => ptr
   endif
end subroutine


! ==============================================================================
! adds a row of real numbers to the table
subroutine add_row_real(t, name, data, format)
   type(table_printer_type), intent(inout) :: t
   character(*), intent(in) :: name    ! name of the row
   real        , intent(in) :: data(:) ! table data
   character(*), intent(in), optional :: format

   character(32) :: fcell
   integer :: i

   call realloc_table(t)

   fcell = DEFAULT_FORMAT; if (present(format)) fcell=format
   t%nr = t%nr + 1
   t%table(0,t%nr) = name
   do i = 1, min(t%nc,size(data))
      write(t%table(i,t%nr),fcell) data(i)
      t%table(i,t%nr) = adjustl(t%table(i,t%nr))
   enddo
end subroutine add_row_real


! ==============================================================================
! adds a row of integer numbers to the table
subroutine add_row_integer(t, name, data, format)
   type(table_printer_type), intent(inout) :: t
   character(*), intent(in) :: name    ! name of the row
   integer     , intent(in) :: data(:) ! table data
   character(*), intent(in), optional :: format

   character(32) :: fcell
   integer :: i

   call realloc_table(t)

   fcell = DEFAULT_FORMAT; if (present(format)) fcell=format
   t%nr = t%nr + 1
   t%table(0,t%nr) = name
   do i = 1, min(t%nc,size(data))
      write(t%table(i,t%nr),fcell) data(i)
      t%table(i,t%nr) = adjustl(t%table(i,t%nr))
   enddo
end subroutine add_row_integer


! ==============================================================================
! adds a row of logical values to the table
subroutine add_row_logical(t, name, data, format)
   type(table_printer_type), intent(inout) :: t
   character(*), intent(in) :: name    ! name of the row
   logical     , intent(in) :: data(:) ! table data
   character(*), intent(in), optional :: format

   character(32) :: fcell
   integer :: i

   call realloc_table(t)

   fcell = DEFAULT_FORMAT; if (present(format)) fcell=format
   t%nr = t%nr + 1
   t%table(0,t%nr) = name
   do i = 1, min(t%nc,size(data))
      write(t%table(i,t%nr),fcell) data(i)
      t%table(i,t%nr) = adjustl(t%table(i,t%nr))
   enddo
end subroutine add_row_logical


! ==============================================================================
! adds a row of character names
subroutine add_row_character(t, name, data)
   type(table_printer_type), intent(inout) :: t
   character(*), intent(in) :: name    ! name of the row
   character(*), intent(in) :: data(:) ! row data

   integer :: i

   call realloc_table(t)
   t%nr = t%nr + 1
   t%table(0,t%nr) = name
   do i = 1, min(t%nc,size(data))
      t%table(i,t%nr) = adjustl(data(i))
   enddo
end subroutine add_row_character

! =============================================================================
! prints the table
subroutine print(t,unit,max_width,cell_width,head_width,transposed)
   type(table_printer_type), intent(in) :: t
   integer, intent(in), optional :: unit       ! unit number for i/o
   integer, intent(in), optional :: max_width  ! maximum width of the cell
   integer, intent(in), optional :: cell_width ! width of each cell
   integer, intent(in), optional :: head_width ! width of the row header
   logical, intent(in), optional :: transposed ! if true, the table is transposed

   integer :: unit_
   integer :: max_width_, head_width_, cell_width_ 
   logical :: transposed_
   integer :: cells_per_line, k

   unit_ = 6 ; if (present(unit))unit_ = unit
   cell_width_ = DEFAULT_CELL_WIDTH  ; if (present(cell_width)) cell_width_ = cell_width
   head_width_ = cell_width_         ; if (present(head_width)) head_width_ = head_width
   max_width_  = DEFAULT_TABLE_WIDTH ; if (present(max_width))  max_width_  = max_width
   transposed_ = .FALSE. ; if (present(transposed)) transposed_ = transposed

   ! calculate max number of cells per line
   cells_per_line = (max_width_-head_width_)/(cell_width_+1)

   if (transposed_) then
     do k = 1, t%nr, cells_per_line
        call print_table_section(transpose(t%table),&
                k,min(t%nr,k+cells_per_line-1),1,t%nc,&
                unit_,cell_width_,head_width_)
     enddo   
   else
     do k = 1, t%nc, cells_per_line
        call print_table_section(t%table,&
                k,min(t%nc,k+cells_per_line-1),1,t%nr,&
                unit_,cell_width_,head_width_)
     enddo
   endif

end subroutine print

! =============================================================================
subroutine print_table_section(table,is,ie,js,je,unit,cell_width,head_width)
   character(*) :: table(0:,0:)
   integer, intent(in) :: is,ie,js,je ! boundaries of the chunk
   integer, intent(in) :: unit       ! unit number for i/o
   integer, intent(in) :: cell_width ! width of each cell
   integer, intent(in) :: head_width ! width of the row header
   
   integer :: i,j,k
   integer :: cells_per_line
   character(128) :: fhead,fcell,ftop,fmid,fbot! format strings
   character(32) :: l, cell
   character :: toprule,midrule,botrule,fieldsep
   
   ! TODO: make all 4 things below optional arguments
   toprule = '='
   midrule = '-'
   botrule = '='
   fieldsep = ' '
   
   ! create formats for the cells, and for the row headers
   write(l,*)cell_width; fcell = '(a'//trim(adjustl(l))//')'
   write(l,*)head_width; fhead = '(a'//trim(adjustl(l))//')'
   
   ! create formats for table separator
   write(l,*)head_width+(cell_width+1)*(ie-is+1)
   l = adjustl(l)
   ftop = "("//trim(l)//"('"//toprule//"'))"
   fmid = "("//trim(l)//"('"//midrule//"'))"
   fbot = "("//trim(l)//"('"//botrule//"'))"

   ! print table header
   write(unit,ftop) ! top-ruler
   write(unit,fhead,advance='no')''
   do i = is,ie
      write(unit,'(a1)',advance='no') fieldsep
      write(unit,fcell,advance='no') table(i,0)
   enddo
   write(unit,*) ! go to the next line
   write(unit,fmid) ! mid-ruler
   ! print table data
   do j = js,je
      write(unit,fhead,advance='no')table(0,j)
      do i = is,ie
         write(unit,'(a1)',advance='no') fieldsep
         if (trim(table(i,j))=='' ) then
            cell = '----'
         else if (len_trim(table(i,j))>cell_width) then
            cell = repeat('*',cell_width)
         else
            cell = table(i,j)
         endif   
         write(unit,fcell,advance='no') cell
      enddo
      write(unit,*) ! go to the next line
   enddo
   write(unit,fbot) ! bottom-ruler
   
end subroutine print_table_section

end module

! ############################################################################
#if 0
program test
   use table_printer_mod
   
   type(table_printer_type) :: t

   call init_with_headers(t,(/'a','b','c','d','e'/))
   call add_row(t,'test1',(/1.0,2.0,3.0,4.0,5.0/))
   call add_row(t,'test2',(/3.5,4.444444,5.6666666e33,1.0,1.0/),format='(e12.4)')
   call add_row(t,'test3',(/4.444444,5.6666666e33,77777.0,2.0,2.0/))
   call add_row(t,'test I',(/11,12,13,14,15/))
   call add_row(t,'test L',(/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE./))
   call add_row(t,'test Char',(/'a1','b2','c3','d4','e5'/))
   call print(t)
   write(*,*)
   call print(t,max_width=70)
   write(*,*)
   call print(t,unit=33,max_width=80,head_width=10)
   call print(t,unit=33,max_width=80,head_width=10,transposed=.TRUE.)
   write(*,*)
   call print(t,head_width=10,cell_width=5)
   write(*,*)
   call print(t,cell_width=10)
   write(*,*)
   call print(t,cell_width=10,transposed=.TRUE.)

   call init_with_headers(t,(/' x','y ','z ','xx','yy','zz'/))
   call add_row(t,'test10',(/1.0,2.0,3.0,4.0,5.0,6.0/))
   call add_row(t,'test11',(/1.0,2.0,3.0,4.0,5.0/))
   call print(t,head_width=10)
   call print(t,head_width=10,transposed=.TRUE.)
   call dealloc(t)
   
end program
#endif
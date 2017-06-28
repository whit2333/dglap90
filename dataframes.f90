! ==============================================================================
! dataframes.f90
!
! by Adam Freese <maxwellsdemon137@gmail.com> copyright 2017
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! This file is part of dglap90.
!
! dglap90 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dglap90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with dglap90.  If not, see <http://www.gnu.org/licenses/>.
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! This file contains the following module
!   - dataframes
!     :: defines a dataframe type, based on the DataFrame type in Python pandas.
!

module dataframes
  ! This module contains the definition of the type dataframe.
  !
  ! dataframe is defined in analogy to DataFrame from Python pandas,
  ! but only has basic functionality.
  ! This type is defiend so that:
  !   1. Matrix columns can be called by a character-valued key.
  !   2. CSV files can be read and written to.
  !
  ! Both of these basic bits of functionality allow smoother interfacing with
  ! scripts that use Python pandas, and this allows for data files to be passed
  ! between programs without having to memorize arbitrary column numbers.
  ! For instance, it doesn't matter whether the column with an up quark
  ! distribution is column 1, or 2, or 8; the up quark column is identified
  ! by the key 'u' in the header.

  use kinds, only : dp

  implicit none
  private

  ! A maximum key length
  integer, parameter :: maxkeylen = 20

  ! A maximum line length for CSV files
  integer, parameter :: maxlinelength = 1024

  ! Silly thing. Can use 'nil' as a failsafe in dataframe%get,
  ! kind of an analogy to using None as the failsafe in Python.
  real(dp), dimension(1), parameter, public :: nil = [ 0.0_dp ]

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! A user-defined type for the data frame
  type, public :: dataframe

    ! A matrix of the actual data. For now, stick with reals
    real(dp), dimension(:,:), allocatable :: matrix(:,:)

    ! The keys for columns
    character(len=maxkeylen), dimension(:), allocatable :: keys

    ! Procedure definitions
    contains
      procedure :: get       => df_get_column
      procedure :: add_col   => df_add_column
      procedure :: read_csv  => io_read_csv
      procedure :: write_csv => io_write_csv
      ! todo :: delete column
      ! todo :: add row(s)
      ! todo :: concatenation

  end type dataframe
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

contains

  ! ============================================================================
  ! ////////////////////////////////////////////////////////////////////////////
  ! "Public" methods (via procedure aliases in dataframe type)
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! ============================================================================

  function df_get_column(this, key, failsafe) result(column) ! alias :: get
    ! Returns column associated with a key.
    ! If the key is not present and an optional failsafe array is provided,
    ! the failsafe array is returned.
    ! If the key is not present and there is no failsafe, the program stops.
    implicit none
    ! I/O
    class(dataframe), intent(in) :: this
    character(len=*), intent(in) :: key
    real(dp), dimension(:), allocatable :: column
    real(dp), dimension(:), intent(in), optional :: failsafe
    ! Intermediates
    integer :: colno
    ! Get the column number
    colno = df_find_key(this, key)
    ! If the key is not found and there's a failsafe, return the failsafe.
    ! Stop the program if there is no failsafe and the search failed.
    if( colno < 0 ) then
      if( present(failsafe) ) then
        column = failsafe
        return
      else
        write(*,*) "Invalid key:", key
        stop
      end if
    end if
    ! Fortran takes care of this allocation automagically
    column = this%matrix(:, colno)
    return
  end function df_get_column

  ! ============================================================================
  ! ============================================================================

  subroutine df_add_column(this, key, column) ! alias :: add_col
    ! Adds column and associates given key with it
    implicit none
    ! Input
    class(dataframe), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), dimension(:), intent(in) :: column
    ! This is either a new allocation, or an append
    if( .not.df_check_alloc(this) ) then
      call df_new_alloc(this, key, column)
    else
      call df_add_alloc(this, key, column)
    end if
    return
  end subroutine df_add_column

  ! ============================================================================
  ! ============================================================================

  subroutine io_read_csv(this, filename) ! alias :: read_csv
    implicit none
    ! I/O
    class(dataframe), intent(inout) :: this
    character(len=*), intent(in) :: filename
    ! Intermediates
    integer :: nrows, ncols
    integer :: ios, iunit
    ! Ensure dataframe is empty before reading into it
    call df_dealloc(this)
    ! Try to open file
    iunit = io_get_free_file_unit()
    open(iunit, file=filename, status='old', iostat=ios)
    if(ios .ne. 0) call exit_failure( "Couldn't open " // trim(filename) )
    ! Get number of rows and columns
    ! First row in file should actually be a list of keys
    nrows = io_get_num_rows(iunit) - 1
    rewind(iunit)
    ncols = io_get_num_cols(iunit)
    rewind(iunit)
    ! First row should be a header, everything else data.
    ! Allocate matrix and keys array appropriately
    allocate( this%matrix(nrows, ncols), this%keys(ncols) )
    ! Copy everything into the dataframe
    call io_read_file( this, iunit )
    ! Close file and be done
    close(iunit)
    return
  end subroutine io_read_csv

  ! ============================================================================
  ! ============================================================================

  subroutine io_write_csv(this, filename) ! alias :: write_csv
    implicit none
    ! I/O
    class(dataframe), intent(in) :: this
    character(len=*), intent(in) :: filename
    ! Intermediates
    integer :: i, nrows
    integer :: ios, iunit
    ! Get the number of rows and columns
    nrows = size(this%matrix, 1)
    ! Try to open file
    iunit = io_get_free_file_unit()
    open(iunit, file=filename, status='replace', iostat=ios)
    if(ios .ne. 0) call exit_failure( "Couldn't open " // trim(filename) )
    ! Write the keys to the top line of the file
    call io_write_key_row(iunit, this%keys)
    ! Write each row of the matrix to file
    do i=1, nrows, 1
      call io_write_real_row( iunit, this%matrix(i,:) )
    end do
    close(iunit)
    return
  end subroutine io_write_csv

  ! ============================================================================
  ! ////////////////////////////////////////////////////////////////////////////
  ! Private (internal) methods
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! ============================================================================

  function df_check_alloc(this) result(isalloc)
    ! Checks whether the dataframe's matrix is allocated yet
    implicit none
    class(dataframe), intent(in) :: this
    logical :: isalloc
    isalloc = .false.
    if( allocated(this%matrix) ) isalloc = .true.
    return
  end function df_check_alloc

  ! ============================================================================
  ! ============================================================================

  subroutine df_new_alloc(this, key, column)
    ! Newly allocates the matrix, adding its first column
    implicit none
    ! Input
    class(dataframe), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), dimension(:), intent(in) :: column
    ! Intermediates
    integer :: nrows, ncols
    ! Don't do new allocation if already allocated
    if( df_check_alloc(this) ) return
    ! Set row and col numbers
    ncols = 1
    nrows = size(column)
    ! Allocate
    allocate( this%matrix(nrows, ncols) )
    allocate( this%keys(ncols) )
    ! Set the first key and the first column
    this%matrix(:,1) = column(:)
    this%keys(1) = key
    ! Be done
    return
  end subroutine df_new_alloc

  ! ============================================================================
  ! ============================================================================

  subroutine df_add_alloc(this, key, column)
    ! Allocates more memory, and appends an additional column
    implicit none
    ! Input
    class(dataframe), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), dimension(:), intent(in) :: column
    ! Intermediates
    real(dp), dimension(:,:), allocatable :: temp_matrix
    character(len=maxkeylen), dimension(:), allocatable :: temp_keys
    integer :: nrows, ncols
    ! Allocate and fill temporary memory
    nrows = size( this%matrix, 1 )
    ncols = size( this%matrix, 2 )
    allocate( temp_matrix(nrows, ncols+1), temp_keys(ncols+1) )
    temp_matrix(:, 1:ncols) = this%matrix(:,1:ncols)
    temp_matrix(:, ncols+1) = column
    temp_keys(1:ncols)      = this%keys(1:ncols)
    temp_keys(ncols+1)      = key
    ! Deallocate the originals
    deallocate( this%matrix, this%keys )
    ! Now move the allocation
    call move_alloc( temp_matrix, this%matrix )
    call move_alloc( temp_keys, this%keys )
    return
  end subroutine df_add_alloc

  ! ============================================================================
  ! ============================================================================

  subroutine df_dealloc(this)
    ! Deallocates dataframe
    implicit none
    ! Input
    class(dataframe), intent(inout) :: this
    ! Deallocate iff allocated
    if( allocated(this%matrix) ) deallocate(this%matrix)
    if( allocated(this%keys) )   deallocate(this%keys)
    return
  end subroutine df_dealloc

  ! ============================================================================
  ! ============================================================================

  function df_find_key(this, key) result(ikey)
    ! Looks for key, returning its index. -1 indicates failure
    implicit none
    ! I/O
    class(dataframe), intent(in) :: this
    character(len=*), intent(in) :: key
    integer :: ikey
    ! Intermediates
    integer :: i
    ! Look for key
    do i=1, size(this%matrix, 2), 1
      if( trim( this%keys(i) ) == trim(key) ) then
        ikey = i
        return
      end if
    end do
    ! Getting this far indicates failure
    ikey = -1
    return
  end function df_find_key

  ! ============================================================================
  ! ////////////////////////////////////////////////////////////////////////////
  ! More private methods, this time related to file I/O
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! ============================================================================

  subroutine io_write_real_row(iunit, row_array)
    implicit none
    ! Input
    integer, intent(in) :: iunit
    real(dp), dimension(:), intent(in) :: row_array
    ! Intermediates
    integer :: N, i
    ! Get the array length and loop over row array elements
    N = size(row_array)
    do i=1, N-1, 1
      write(iunit, 666, advance='no') row_array(i), ','
    end do
    ! No column after the last element, but newline instead
    write(iunit,*) row_array(N)
    ! Be done
    return
666 format(ES20.10,A1)
  end subroutine io_write_real_row

  ! ============================================================================
  ! ============================================================================

  subroutine io_write_key_row(iunit, row_array)
    implicit none
    ! Input
    integer, intent(in) :: iunit
    character(len=*), dimension(:), intent(in) :: row_array
    ! Intermediates
    integer :: N, i
    ! Get the array length and loop over row array elements
    N = size(row_array)
    do i=1, N-1, 1
      write(iunit, 666, advance='no') trim(adjustl( row_array(i) )), ','
    end do
    ! No column after the last element, but newline instead
    write(iunit,*) row_array(N)
    ! Be done
    return
666 format(A,A1)
  end subroutine io_write_key_row

  ! ============================================================================
  ! ============================================================================

  subroutine io_read_file(this, iunit)
    implicit none
    ! Input
    class(dataframe), intent(inout) :: this
    integer, intent(in) :: iunit
    ! Intermediates
    character(len=maxlinelength) :: line
    integer :: irow, ios
    integer :: icnt ! debug
    ! Loop
    irow = 0
    icnt = 0
    do
      icnt = icnt + 1
      read(iunit, '(A)', iostat=ios) line
      if( .not.io_is_comment(line) ) then
        backspace(iunit)
        ! Read into keys or matrix, depending on irow
        if(irow==0) read(iunit,*) this%keys
        if(irow >0) read(iunit,*) this%matrix(irow,:)
        irow = irow + 1
      end if
      ! Quit when we've exceeded the number of rows
      if( irow > size(this%matrix,1) ) return
    end do
  end subroutine io_read_file

  ! ============================================================================
  ! ============================================================================

  function io_get_num_rows(iunit) result(numrows)
    ! Obtains the number of columns in a CSV file
    implicit none
    ! I/O
    integer, intent(in) :: iunit
    integer :: numrows
    ! Intermediates
    character(len=maxlinelength) :: line
    integer :: ios
    ! Keep reading until end of file
    numrows = 0
    do
      read(iunit, '(A)', iostat=ios) line
      ! Leave once there are no more lines to read
      if( ios .ne. 0 ) return
      ! Check whether line is a comment
      if( .not.io_is_comment(line) ) numrows = numrows + 1
    end do
  end function io_get_num_rows

  ! ============================================================================
  ! ============================================================================

  function io_get_num_cols(iunit) result(numcols)
    ! Obtains the number of columns in a CSV file
    implicit none
    ! I/O
    integer, intent(in) :: iunit
    integer :: numcols
    ! Intermediates
    character(len=maxlinelength) :: line
    integer :: ios
    ! Keep reading until we find an acceptable line
    do
      read(iunit, '(A)', iostat=ios) line
      ! Quit if we don't find anything.
      if( ios .ne. 0 ) call exit_failure("Couldn't parse file.")
      ! Check whether line is a comment
      if( .not.io_is_comment(line) ) then
        ! For a CSV file, column count is commas+1
        ! Unless user put comma in a string, which would be silly...
        numcols = io_count_commas(line) + 1
        return
      end if
    end do
  end function io_get_num_cols

  ! ============================================================================
  ! ============================================================================

  function io_is_comment(line) result(tof)
    ! Checks if given line from a file is a comment
    implicit none
    ! I/O
    character(len=*), intent(in) :: line
    logical :: tof
    ! Intermediates
    character(len=len(line)) :: eff_line
    ! Innocent until proven guilty
    tof = .false.
    ! Get rid of trailing whitespace on both sides
    eff_line = trim( adjustl(line) )
    ! If the line is empty or just whitespace, count it as a comment
    if( len_trim(eff_line) == 0 ) then
      tof = .true.
      return
    end if
    ! If line starts with a hash, it's a comment
    if( eff_line(1:1) == '#' ) tof = .true.
    return
  end function io_is_comment

  ! ============================================================================
  ! ============================================================================

  function io_count_commas(line) result(commacount)
    ! Counts number of commas in a string
    implicit none
    character(len=*), intent(in) :: line
    integer :: commacount
    ! Intermediates
    integer :: i
    ! Start at zero
    commacount = 0
    do i=1, len(line), 1
      if( line(i:i) == ',' ) commacount = commacount + 1
    end do
    return
  end function io_count_commas

  ! ============================================================================
  ! ============================================================================

  function io_get_free_file_unit() result(iunit)
    implicit none
    ! I/O
    integer :: iunit
    ! Intermediates
    integer ios
    logical :: lopen
    ! Look for a free file unit
    do iunit = 1, 99, 1
       ! Skip 5 and 6
       if( iunit.ne.5 .and. iunit.ne.6 ) then
         ! Return with the iunit if status is OK and there are no errors.
          inquire( unit = iunit, opened = lopen, iostat = ios )
          if ( ios.eq.0 .and. .not.lopen ) return
       end if
    end do
    ! If we got this far, failure; the program dies.
    call exit_failure("Could not find a free file unit.")
  end function io_get_free_file_unit

  ! ============================================================================
  ! ////////////////////////////////////////////////////////////////////////////
  ! And more private methods
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! ============================================================================

  subroutine exit_failure(exit_message)
    implicit none
    character(len=*), intent(in) :: exit_message
    write(*,*) "FATAL: " // exit_message
    stop
  end subroutine exit_failure

end module dataframes

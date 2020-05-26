module input_read_module
#include "petsc/finclude/petscsys.h"
  use petscsys
!
  implicit none
!
  public :: StringToUpper, &
            StringToLower, &
            InputReadString, &
            InputReadWord, &
            StringAdjustl, &
            InputCreate
  integer, parameter, public :: MAXSTRINGLENGTH = 512
  PetscInt, parameter, public :: MAXWORDLENGTH = 32

contains
!
subroutine StringAdjustl(string)
  ! 
  ! Left adjusts a string by removing leading spaces and tabs.
  ! This subroutine is needed because the adjustl() Fortran 90
  ! intrinsic will not remove leading tabs.
  ! 
  ! Author: Richard Tran Mills
  ! Date: 9/21/2010
  ! 

  implicit none

  character(len=*) :: string

  PetscInt :: i
  PetscInt :: string_length
  character(len=1), parameter :: tab = achar(9)

  ! We have to manually convert any leading tabs into spaces, as the 
  ! adjustl() intrinsic does not eliminate leading tabs.
  i=1
  string_length = len_trim(string)
  do while((string(i:i) == ' ' .or. string(i:i) == tab) .and. &
           i <= string_length)
    if (string(i:i) == tab) string(i:i) = ' '
    i=i+1
  enddo
  ! adjustl() will do what we want, now that tabs are removed.
  string = adjustl(string)

end subroutine StringAdjustl


subroutine StringToUpper(string)
  ! 
  ! converts lowercase characters in a card to uppercase
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
      string(i:i) = achar(iachar(string(i:i)) - 32)
    endif
  enddo

end subroutine StringToUpper

! ************************************************************************** !

subroutine StringToLower(string)
  ! 
  ! converts uppercase characters in a card to lowercase
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
      string(i:i) = achar(iachar(string(i:i)) + 32)
    endif
  enddo

end subroutine StringToLower
! ************************************************************************** !


subroutine InputReadString(fid,buf,ierr)
  ! 
  ! Reads a string (strlen characters long) from a
  ! file while avoiding commented or skipped lines.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 


  implicit none

  integer :: fid
  character(len=MAXSTRINGLENGTH) ::  tempstring
  character(len=MAXWORDLENGTH) :: word
  character(len=*)             :: buf
  integer :: i
  integer :: skip_count
  integer :: ierr
  ierr = 0

  word = ''

  do
    read(fid,'(a512)',iostat=ierr) buf
    call StringAdjustl(buf)

    if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle

    tempstring = buf
    call InputReadWord(tempstring,word,.true.,ierr)
    call StringToUpper(word)

    if (word(1:4) == 'SKIP') then
      ! to avoid keywords that start with SKIP 
      if (len_trim(word) > 4) then
        exit
      endif
      skip_count = 1
      do
        read(fid,'(a512)',iostat=ierr) tempstring

        call InputReadWord(tempstring,word,PETSC_FALSE,ierr)
        call StringToUpper(word)
        if (word(1:4) == 'SKIP') skip_count = skip_count + 1
        if (word(1:4) == 'NOSK') then
          skip_count = skip_count - 1
          if (skip_count == 0) exit
        endif
      enddo
      if (ierr /= 0) exit
    else if (word(1:1) /= ' ' .and. word(1:4) /= 'NOSK') then
      exit
    endif
  enddo

  ! Check for comment midway along a string
  if (.not.InputError(ierr)) then
    tempstring = buf
    buf = repeat(' ',MAXSTRINGLENGTH)
    do i=1,len_trim(tempstring)
      if (tempstring(i:i) /= '#' .and. tempstring(i:i) /= '!') then
        buf(i:i) = tempstring(i:i)
      else
        exit
      endif
    enddo
  endif

end subroutine InputReadString

  
subroutine InputReadWord(string, word, return_blank_error, ierr)
  ! 
  ! reads and removes a word (consecutive characters) from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  character(len=*) :: string
  character(len=*) :: word
  PetscBool :: return_blank_error
  PetscErrorCode :: ierr

  PetscInt :: i, begins, ends, length
  character(len=1), parameter :: tab = achar(9), backslash = achar(92)

  if (ierr /= 0) return
  word = ''

  length = len_trim(string)

  if (length == 0) then
    if (return_blank_error) then
      ierr = 1
    else
      ierr = 0
    endif
    return
  else
    ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while((string(i:i) == ' ' .or. string(i:i) == ',' .or. &
             string(i:i) == tab) .and. i <= length)
      i=i+1
    enddo

    if (i > length) then
      if (return_blank_error) then
        ierr = 1
      else
        ierr = 0
      endif
      return
    endif

    begins=i
    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab .and. &
              (i == begins .or. string(i:i) /= backslash))
      i=i+1
    enddo

    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

    ! Copy (ends-begins) characters to 'word'
    word = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif
end subroutine InputReadWord
!
subroutine InputCreate(fid,path,filename)
  ! 
  ! Allocates and initializes a new Input object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 


  implicit none

  integer :: fid
  character(len=*) :: path
  character(len=*) :: filename
  character(len=MAXSTRINGLENGTH) :: tfilename
  character(len=MAXSTRINGLENGTH)             :: buf
  character(len=MAXSTRINGLENGTH)             :: err_buf
  character(len=MAXSTRINGLENGTH)             :: err_buf2

  integer :: istatus
  integer :: islash
  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: local_path
  character(len=MAXSTRINGLENGTH) :: full_path
  logical, parameter :: back = .true.

  fid = fid
  path = ''
!  filename = ''
  tfilename = ''
  ierr = 0
  buf = ''
  err_buf = ''
  err_buf2 = ''

  ! split the filename into a path and filename
                              ! backwards search
  islash = index(filename,'/',back)
  if (islash > 0) then
    path(1:islash) = filename(1:islash)
    tfilename(1:len_trim(filename)-islash) = &
      filename(islash+1:len_trim(filename))
  else
    tfilename = filename
  endif


  full_path = trim(path) // trim(tfilename)
  open(unit=fid,file=full_path,status="old",iostat=istatus)
!  if (istatus /= 0) then
!    if (len_trim(full_path) == 0) full_path = '<blank>'
!    io_buffer = 'File: "' // trim(full_path) // '" not found.'
!  endif
end subroutine InputCreate

subroutine InputReadInt(buf,int,ierr)
  ! 
  ! reads and removes an integer value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: int

  character(len=*)             :: buf
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found
  integer :: ierr

  ierr = 0 
  found = .false.

  if (.not.found) then
    call InputReadWord(buf,word,.true.,ierr)

    if (.not.InputError(ierr)) then
      read(word,*,iostat=ierr) int
    endif
  endif

end subroutine InputReadInt

subroutine InputReadDouble(buf, double,ierr)
  ! 
  ! reads and removes a real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscReal :: double
  character(len=*)             :: buf
 
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found
  integer :: ierr
  ierr = 0
  found = .false.
  if(.not. found) then
    call InputReadWord(buf,word,.true.,ierr)

    if (.not.InputError(ierr)) then
        read(word,*,iostat=ierr) double
    endif
  endif

end subroutine InputReadDouble

logical function StringCompare(string1,string2,n)
  ! 
  ! compares two strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: i, n
  character(len=n) :: string1, string2

  do i=1,n
    if (string1(i:i) /= string2(i:i)) then
      StringCompare = .false.
      return
    endif
  enddo

  StringCompare = .true.
  return

end function StringCompare

subroutine InputReadNChars(string, chars, n, return_blank_error, ierr)
  ! 
  ! reads and removes a specified number of characters from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/00
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.

  PetscInt :: i, n, begins, ends
  character(len=n) :: chars
  integer :: ierr
  character(len=1), parameter :: tab = achar(9), backslash = achar(92)

  if (inputError(ierr)) return
  ! Initialize character string to blank.
  chars(1:n) = repeat(' ',n)

  ierr = len_trim(string)
  if (.not.inputError(ierr)) then
    if (return_blank_error) then
      ierr = 1
    else
      ierr = 0
    endif
    return
  else
    ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while(string(i:i) == ' ' .or. string(i:i) == tab)
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab  .and. &
              (i == begins .or. string(i:i) /= backslash))
      i=i+1
    enddo

    ends=i-1

    if (ends-begins+1 > n) then ! string read is too large for 'chars'
      ierr = 1
      return
    endif

    ! Copy (ends-begins) characters to 'chars'
    chars = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine InputReadNChars

function InputError(ierr)
  !
  ! Returns true if an error has occurred
  !
  ! Author: Glenn Hammond
  ! Date: 12/10/08
  !

  implicit none

  PetscErrorCode :: ierr

  PetscBool :: InputError

  if (ierr == 0) then
    InputError = PETSC_FALSE
  else
    InputError = PETSC_TRUE
  endif

end function InputError
    
end module input_read_module

!
! Copyright (C) 2001-2004 Carlo Cavazzoni and PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION matches( string1, string2 )  
  !-----------------------------------------------------------------------
  !
  ! ... .TRUE. if string1 is contained in string2, .FALSE. otherwise
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*), INTENT(IN) :: string1, string2
  LOGICAL                       :: matches
  INTEGER                       :: len1, len2, l  
  !
  !
  len1 = LEN_TRIM( string1 )  
  len2 = LEN_TRIM( string2 )  
  !
  DO l = 1, ( len2 - len1 + 1 )  
     !   
     IF ( string1(1:len1) == string2(l:(l+len1-1)) ) THEN  
        !
        matches = .TRUE.  
        !
        RETURN  
        !
     END IF
     !
  END DO
  !
  matches = .FALSE.
  ! 
  RETURN
  !
END FUNCTION matches
!
!-----------------------------------------------------------------------
FUNCTION imatches( string1, string2 )
  !-----------------------------------------------------------------------
  !
  ! ... .TRUE. if string1 is contained in string2, .FALSE. otherwise
  !   *** case insensitive ***
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*), INTENT(IN) :: string1, string2
  CHARACTER(LEN=len(string1))   :: aux1
  CHARACTER(LEN=len(string2))   :: aux2
  LOGICAL                       :: imatches
  LOGICAL, EXTERNAL             :: matches
  !
  aux1 = string1
  aux2 = string2
  ! 
  CALL up2lw(aux1)
  CALL up2lw(aux2)
  !
  imatches = matches(aux1, aux2)
  !
  RETURN
  !
END FUNCTION imatches
!
!-----------------------------------------------------------------------
SUBROUTINE up2lw(word)
  !-----------------------------------------------------------------------
  ! convert a word to lower case 
  character (len=*) , intent(in out) :: word 
  integer :: i,ic,nlen 
  nlen = len(word) 
  do i=1,nlen 
     ic = ichar(word(i:i)) 
     if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32) 
  end do
END SUBROUTINE up2lw
!
!-----------------------------------------------------------------------
SUBROUTINE lw2up(word)
  !-----------------------------------------------------------------------
  ! convert a word to upper case 
  character (len=*) , intent(in out) :: word 
  integer :: i,ic,nlen 
  nlen = len(word) 
  do i=1,nlen 
     ic = ichar(word(i:i)) 
     if (ic >= 97 .and. ic < 122) word(i:i) = char(ic-32) 
  end do
END SUBROUTINE lw2up




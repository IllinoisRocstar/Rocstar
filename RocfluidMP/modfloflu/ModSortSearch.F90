! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
! ******************************************************************************
!
! Purpose: Collection of procedures to sort and search arrays.
!
! Description: none
!
! Notes:
!   1. The procedures in this modules have been collected from various sources. 
!      See remarks in the respective procedures.
!
! ******************************************************************************
!
! $Id: ModSortSearch.F90,v 1.18 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModSortSearch

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Local variables
! ==============================================================================
 
  INTEGER, PARAMETER :: ELEMENT_NOT_FOUND = -1 ! must be < 0 
  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: & 
    RCSIdentString = '$RCSfile: ModSortSearch.F90,v $ $Revision: 1.18 $'

! ******************************************************************************
! Module subroutines
! ******************************************************************************

  CONTAINS
  
  
  
  
! ******************************************************************************
!
! Procedure: BinarySearchInteger.F90
!
! Purpose: Search array with elements in ascending order using binary search.
!
! Description: See appropriate textbooks.
!
! Input:
!   a   Array with elements in ascending order
!   n   Size of array
!   v   Value which is being searched for
! 
! Output:
!   i   Location of v in a
!   j   Lower index if value was not found (optional)
!
! Notes: 
!   1. If v is not in a, the error flag ELEMENT_NOT_FOUND is returned.
!   2. The optional argument j gives lower index if value was not found, and 
!      CRAZY_VALUE_INT if the value was found. It is important to note that 
!      this lower index can exceed n, the dimension of the array being searched.
!      For this reason, it is important to make sure j is checked before it is
!      used to access an array also dimensioned with n.
!
! ******************************************************************************

SUBROUTINE BinarySearchInteger(a,n,v,i,j)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Parameters
! ==============================================================================   
  
  INTEGER, INTENT(IN) :: n,v
  INTEGER, INTENT(IN) :: a(n)
  INTEGER, INTENT(OUT) :: i
  INTEGER, INTENT(INOUT), OPTIONAL :: j
  
! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: il,im,iu
  
! ******************************************************************************
! Start
! ******************************************************************************

  il = 1 ! initialise lower limit
  iu = n ! initialise upper limit

  DO 
    im = (iu + il)/2 ! compute midpoint
    
    IF ( v < a(im) ) THEN ! eliminate upper half
      iu = im - 1 
    ELSE IF ( v > a(im) ) THEN ! eliminate lower half
      il = im + 1 
    ELSE ! found value
      i = im
      
      IF ( PRESENT(j) .EQV. .TRUE. ) THEN 
        j = CRAZY_VALUE_INT
      END IF ! PRESENT      
      
      EXIT    
    END IF ! v
    
    IF ( iu < il ) THEN ! element not found
      i = ELEMENT_NOT_FOUND
      
      IF ( PRESENT(j) .EQV. .TRUE. ) THEN 
        j = il
      END IF ! PRESENT 
      
      EXIT
    END IF ! iu
  END DO ! <empty> 

! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE BinarySearchInteger
  
  
  
  


! ******************************************************************************
!
! Purpose: Cycle list of integers until condition of equality is met.
!
! Description: None.
!
! Input:
!   a   Array containing n elements
!   na  Number of elements in array a
!   i   Location at a where value must be located
!   v   Target value in a at location i
!
! Output:
!   a   Array cycled such that a(i) = v
!
! Notes: 
!   1. This routine does not check that v actually exists in a, so it is up to
!      you to make sure before calling this routine!
!
! ******************************************************************************

  SUBROUTINE CycleList(a,na,i,v)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   

    INTEGER, INTENT(IN) :: i,na,v
    INTEGER, INTENT(INOUT) :: a(na)

! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: j,k,w

! ******************************************************************************
!   Start
! ******************************************************************************  

    DO j = 1,na-1
      IF ( a(i) /= v ) THEN 
        w = a(na)

        DO k = na,2,-1
          a(k) = a(k-1)  
        END DO ! k

        a(1) = w
      ELSE 
        EXIT         
      END IF ! a            
    END DO ! j

  END SUBROUTINE CycleList
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Find duplicated elements of sorted integer arrays a and b and 
!   store in array c, while keeping all elements in original arrays.
!
! Description: None.
!
! Input:
!   a           Array 
!   na          Number of elements in array a
!   b           Array, must be sorted in increasing order
!   nb          Number of elements in array b
!   c           Array 
!   ncMax       Size of array c
!
! Output:
!   c           Array c with duplicated elements on a and b
!   nc          Number of elements in array c
!   errorFlag   Error flag
!
! Notes: 
!   1. List b must be sorted in increasing order on entry. List a, strictly
!      speaking must not be sorted. If a is not sorted, then the array c will
!      not be sorted either.
!
! ******************************************************************************

  SUBROUTINE FindCommonSortedIntegers(a,na,b,nb,c,ncMax,nc,errorFlag)
  
    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   
    
    INTEGER, INTENT(IN) :: na,nb,ncMax    
    INTEGER, INTENT(OUT) :: nc,errorFlag
    INTEGER, INTENT(IN) :: a(na),b(nb)
    INTEGER, INTENT(OUT) :: c(ncMax) 
    
! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: ia,iLoc

! ******************************************************************************
!   Start
! ******************************************************************************  

    errorFlag = ERR_NONE

    nc = 0 

    ia = 1

    emptyLoop: DO
      IF ( nb > 1 ) THEN  
        CALL BinarySearchInteger(b(1:nb),nb,a(ia),iLoc)
      ELSE 
        IF ( b(1) == a(ia) ) THEN 
          iLoc = ELEMENT_NOT_FOUND + 1 ! Set to anything but ELEMENT_NOT_FOUND
        ELSE 
          iLoc = ELEMENT_NOT_FOUND
        END IF ! b(1)
      END IF ! nb

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN       
        nc = nc + 1
                
        IF ( nc > ncMax ) THEN 
          errorFlag = ERR_NONE - 1 ! Set to anything but ERR_NONE
          EXIT emptyLoop
        END IF ! ic
        
        c(nc) = a(ia) 
      END IF ! iLoc

      ia = ia + 1
      
      IF ( ia > na ) THEN 
        EXIT emptyLoop
      END IF ! ia
    END DO emptyLoop
  
  END SUBROUTINE FindCommonSortedIntegers  
  
  
  
  


! ******************************************************************************
!
! Purpose: Merge elements of two sorted elements of integer arrays a of size n.
!
! Description: None.
!
! Input:
!   na  Number of elements in array a
!   nb  Number of elements in array b 
!   a   Array containing sorted elements
!   b   Array containing sorted elements
!   nm  Number of elements in array m
!   m   Array (empty)
!
! Output:
!   im  Number of elements in array m
!   m   Array containing merged elements
!
! Notes: 
!   1. This procedure results in a unique list in the sense that each element
!      in the merged list only appears once. 
!   2. It is easier to remove duplicate elements from the lists before merging 
!      them than dealing with duplicate elements during merging.
!
! ******************************************************************************

  SUBROUTINE MergeSortedIntegers(global,na,nb,a,b,nm,im,m)

    IMPLICIT NONE
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   
    
    INTEGER, INTENT(IN) :: na,nb,nm
    INTEGER, INTENT(IN) :: a(na),b(nb)
    TYPE(t_global), POINTER :: global
    
    INTEGER, INTENT(OUT) :: im
    INTEGER, INTENT(OUT) :: m(nm)

! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: i,ia,ib,na2,nb2,tbi
    INTEGER :: a2(na),b2(nb)

! ******************************************************************************
!   Start
! ******************************************************************************         
                    
! Eliminate duplicate elements from lists

    a2(1:na) = a(1:na)
    b2(1:nb) = b(1:nb)

    CALL SimplifySortedIntegers(a2,na,na2)
    CALL SimplifySortedIntegers(b2,nb,nb2)   
                      
! Initialize
        
    m(1:nm) = 0
    
    ia = 1
    ib = 1
    im = 0
    
! Traverse both sorted lists and insert smaller value into m    
    
    DO 
      IF ( ia <= na2 .AND. ib <= nb2 ) THEN 
        IF ( a2(ia) > b2(ib) ) THEN
          tbi = b2(ib) 
          ib  = ib + 1
        ELSE 
          tbi = a2(ia) 
          ia  = ia + 1
        END IF ! a2(ia)
        
        IF ( im == 0 ) THEN 
          im = 1
          m(im) = tbi            
        ELSE 
          IF ( im < nm ) THEN 
            IF ( tbi /= m(im) ) THEN 
              im = im + 1              
              m(im) = tbi
            END IF ! tbi          
          ELSE
            CALL ErrorStop(global,ERR_MERGE_SORTED,__LINE__)
          END IF ! im        
        END IF ! im        
      ELSE 
        EXIT
      END IF ! ia ...
    
    END DO ! <empty>

! Now that one of the original lists was traversed completely, fill in 
! the remainder of m from other sorted list

    IF ( ia <= na2 ) THEN 
      DO i = ia,na2
        IF ( a2(i) /= m(im) ) THEN ! do not duplicate element
          im = im + 1
          m(im) = a2(i)
        END IF ! a2(i)      
      END DO ! i
    ELSE
      DO i = ib,nb2
        IF ( b2(i) /= m(im) ) THEN ! do not duplicate element
          im = im + 1
          m(im) = b2(i)
        END IF ! b2(i)
      END DO ! i
    END IF ! ia

  END SUBROUTINE MergeSortedIntegers

  
  
 
  
  
! ******************************************************************************
!
! Procedure: QuickSortInteger.F90
!
! Purpose: Sort elements of integer array a of size n into ascending order.
!
! Description: Uses quicksort method. 
!
! Input:
!   a   Array containing unsorted elements
!   n   Number of elements in array a
! 
! Output:
!   a   Array containing sorted elements
!
! Notes:
!   1. Taken from WWW, cannot remember where... Seems to originate from 
!      Nicklaus Wirths book, see remark below.
!   2. No modifications to original, apart from a few cosmetic changes and 
!      from deletion of second array which was also being sorted along with a.
!
! ******************************************************************************

  SUBROUTINE QuickSortInteger(a,n)

! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTHS PASCAL
! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(INOUT) :: a(n)

! Local Variables

    INTEGER :: i, j, k, l, r, s, stackl(50), stackr(50)
    INTEGER :: w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

 10 CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

 20 CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

! REPEAT UNTIL I > J.

    DO
      DO
        IF (a(i).LT.x) THEN ! Search from lower end
          i = i + 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      DO
        IF (x.LT.a(j)) THEN ! Search from upper end
          j = j - 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      IF (i.LE.j) THEN ! Swap positions i & j
        w = a(i)
        a(i) = a(j)
        a(j) = w
        i = i + 1
        j = j - 1
        IF (i.GT.j) EXIT
      ELSE
        EXIT
      END IF
    END DO

    IF (j-l.GE.r-i) THEN
      IF (l.LT.j) THEN
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      END IF
      l = i
    ELSE
      IF (i.LT.r) THEN
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      END IF
      r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN

  END SUBROUTINE QuickSortInteger





! ******************************************************************************
!
! Procedure: QuickSortIntegerInteger.F90
!
! Purpose: Sort elements of integer array a of size n into ascending order and 
!   sort b along with it.
!
! Description: Uses quicksort method. 
!
! Input:
!   a   Array containing unsorted elements
!   b   Integer array
!   n   Number of elements in array a
! 
! Output:
!   a   Array containing sorted elements
!   b   Integer array
!
! Notes:
!   1. Taken from WWW, cannot remember where... Seems to originate from 
!      Nicklaus Wirths book, see remark below.
!   2. No modifications to original, apart from a few cosmetic changes and 
!      from deletion of second array which was also being sorted along with a.
!
! ******************************************************************************

  SUBROUTINE QuickSortIntegerInteger(a,b,n)

! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTHS PASCAL
! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(INOUT) :: a(n),b(n)

! Local Variables

    INTEGER :: i, j, k, l, r, s, stackl(50), stackr(50)
    INTEGER :: v, w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

 10 CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

 20 CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

! REPEAT UNTIL I > J.

    DO
      DO
        IF (a(i).LT.x) THEN ! Search from lower end
          i = i + 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      DO
        IF (x.LT.a(j)) THEN ! Search from upper end
          j = j - 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO
      
      IF (i.LE.j) THEN ! Swap positions i & j
        w = a(i)
        a(i) = a(j)
        a(j) = w
        v = b(i)
        b(i) = b(j)
        b(j) = v
        i = i + 1
        j = j - 1
        IF (i.GT.j) EXIT
      ELSE
        EXIT
      END IF
    END DO

    IF (j-l.GE.r-i) THEN
      IF (l.LT.j) THEN
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      END IF
      l = i
    ELSE
      IF (i.LT.r) THEN
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      END IF
      r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN

  END SUBROUTINE QuickSortIntegerInteger






! ******************************************************************************
!
! Purpose: Sort elements of double-precision array a of size n into ascending 
!   order.
!
! Description: Uses quicksort method. 
!
! Input:
!   a   Array containing unsorted elements
!   n   Number of elements in array a
! 
! Output:
!   a   Array containing sorted elements
!
! Notes:
!   1. Taken from WWW, cannot remember where... Seems to originate from 
!      Nicklaus Wirths book, see remark below.
!   2. No modifications to original, apart from a few cosmetic changes and 
!      from deletion of second array which was also being sorted along with a.
!
! ******************************************************************************

  SUBROUTINE QuickSortRFREAL(a,n)

! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTHS PASCAL
! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(KIND=RFREAL), INTENT(INOUT) :: a(n)

! Local Variables

    INTEGER :: i, j, k, l, r, s, stackl(50), stackr(50)
    REAL(KIND=RFREAL) :: w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

 10 CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

 20 CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

! REPEAT UNTIL I > J.

    DO
      DO
        IF (a(i).LT.x) THEN ! Search from lower end
          i = i + 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      DO
        IF (x.LT.a(j)) THEN ! Search from upper end
          j = j - 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      IF (i.LE.j) THEN ! Swap positions i & j
        w = a(i)
        a(i) = a(j)
        a(j) = w
        i = i + 1
        j = j - 1
        IF (i.GT.j) EXIT
      ELSE
        EXIT
      END IF
    END DO

    IF (j-l.GE.r-i) THEN
      IF (l.LT.j) THEN
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      END IF
      l = i
    ELSE
      IF (i.LT.r) THEN
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      END IF
      r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN

  END SUBROUTINE QuickSortRFREAL





! ******************************************************************************
!
! Purpose: Sort elements of double-precision array a of size n into ascending 
!   order and sort b along with it.
!
! Description: Uses quicksort method. 
!
! Input:
!   a   Array containing unsorted elements
!   b   Integer array
!   n   Number of elements in array a
! 
! Output:
!   a   Array containing sorted elements
!   b   Integer array containing elements in same sorted order as a
!
! Notes:
!   1. Taken from WWW, cannot remember where... Seems to originate from 
!      Nicklaus Wirths book, see remark below.
!   2. No modifications to original, apart from a few cosmetic changes and 
!      from deletion of second array which was also being sorted along with a.
!
! ******************************************************************************

  SUBROUTINE QuickSortRFREALInteger(a,b,n)

! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTHS PASCAL
! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(INOUT) :: b(n)    
    REAL(KIND=RFREAL), INTENT(INOUT) :: a(n)

! Local Variables

    INTEGER :: i, j, k, l, r, s, stackl(50), stackr(50), v
    REAL(KIND=RFREAL) :: w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

 10 CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

 20 CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

! REPEAT UNTIL I > J.

    DO
      DO
        IF (a(i).LT.x) THEN ! Search from lower end
          i = i + 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      DO
        IF (x.LT.a(j)) THEN ! Search from upper end
          j = j - 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      IF (i.LE.j) THEN ! Swap positions i & j
        w = a(i)
        a(i) = a(j)
        a(j) = w
        v = b(i)
        b(i) = b(j)
        b(j) = v 
        i = i + 1
        j = j - 1
        IF (i.GT.j) EXIT
      ELSE
        EXIT
      END IF
    END DO

    IF (j-l.GE.r-i) THEN
      IF (l.LT.j) THEN
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      END IF
      l = i
    ELSE
      IF (i.LT.r) THEN
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      END IF
      r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN

  END SUBROUTINE QuickSortRFREALInteger







! ******************************************************************************
!
! Purpose: Remove duplicated elements of sorted integer arrays a of size na.
!
! Description: None.
!
! Input:
!   a   Array 
!   na  Number of elements in array a
!   b   Array, must be sorted in increasing order
!   nb  Number of elements in array b
!
! Output:
!   b   Array without elements shared with a
!   nb2 Number of elements in array b which are not shared with a
!
! Notes: 
!   1. List b must be sorted in increasing order on entry. 
!
! ******************************************************************************

  SUBROUTINE RemoveCommonSortedIntegers(a,na,b,nb,nb2)
  
    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   
    
    INTEGER, INTENT(IN) :: na
    INTEGER, INTENT(IN) :: nb
    INTEGER, INTENT(OUT) :: nb2
    INTEGER, INTENT(IN) :: a(na)       
    INTEGER, INTENT(INOUT) :: b(nb) 
    
! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: ia,ib2,iLoc

! ******************************************************************************
!   Start
! ******************************************************************************  

    nb2 = nb

    DO ia = 1,na
      CALL BinarySearchInteger(b(1:nb2),nb2,a(ia),iLoc)
      
      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN
        CALL RemoveInteger(b,nb2,iLoc)
      END IF ! iLoc
    END DO ! ia
  
  END SUBROUTINE RemoveCommonSortedIntegers






! ******************************************************************************
!
! Purpose: Remove duplicated elements of sorted integer arrays a and b and 
!   store in array c, while keeping non-duplicated elements in original arrays.
!
! Description: None.
!
! Input:
!   a           Array 
!   na          Number of elements in array a
!   b           Array, must be sorted in increasing order
!   nb          Number of elements in array b
!   c           Array 
!   ncMax       Size of array c
!
! Output:
!   a           Array a without duplicated elements
!   na2         Number of elements in array a
!   b           Array b without duplicated elements
!   nb2         Number of elements in array b
!   c           Array c with duplicated elements on a and b
!   nc          Number of elements in array c
!   errorFlag   Error flag indicating whether routine completed successfully.
!               Any non-zero value indicates failure.
!
! Notes: 
!   1. List b must be sorted in increasing order on entry. List a, strictly
!      speaking must not be sorted. If a is not sorted, then the array c will
!      not be sorted either.
!
! ******************************************************************************

  SUBROUTINE RemoveCommonSortedIntegersFancy(a,na,na2,b,nb,nb2,c,ncMax,nc, & 
                                             errorFlag)
  
    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   
    
    INTEGER, INTENT(IN) :: na,nb,ncMax    
    INTEGER, INTENT(OUT) :: errorFlag,na2,nb2,nc
    INTEGER, INTENT(INOUT) :: a(na),b(nb)
    INTEGER, INTENT(OUT) :: c(ncMax) 
    
! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: ia,ib2,iLoc

! ******************************************************************************
!   Start
! ******************************************************************************  

    errorFlag = ERR_NONE

    na2 = na
    nb2 = nb
    nc  = 0 

    ia = 1

    emptyLoop: DO
      IF ( nb2 > 1 ) THEN  
        CALL BinarySearchInteger(b(1:nb2),nb2,a(ia),iLoc)
      ELSE 
        IF ( b(1) == a(ia) ) THEN 
          iLoc = ELEMENT_NOT_FOUND + 1 ! Set to anything but ELEMENT_NOT_FOUND
        ELSE 
          iLoc = ELEMENT_NOT_FOUND
        END IF ! b(1)
      END IF ! nb

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN       
        nc = nc + 1
        
        IF ( nc > ncMax ) THEN 
          errorFlag = ERR_NONE - 1 ! Set to anything but ERR_NONE
          EXIT emptyLoop
        END IF ! ic
        
        c(nc) = a(ia) 
        
        CALL RemoveInteger(a(1:na2),na2,ia)
        CALL RemoveInteger(b(1:nb2),nb2,iLoc)
      ELSE 
        ia = ia + 1
      END IF ! iLoc
      
      IF ( ia > na2 ) THEN 
        EXIT emptyLoop
      END IF ! ia
    END DO emptyLoop
  
  END SUBROUTINE RemoveCommonSortedIntegersFancy





! ******************************************************************************
!
! Purpose: Remove entry from integer array and set last element to crazy value.
!
! Description: None.
!
! Input:
!   a    Array 
!   na   Number of elements in array a
!   iLoc Position of entry in array a to be removed 
!
! Output:
!   a    Array with iLoc-th element removed
!   na   Number of elements in array a after removing iLoc-th element
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RemoveInteger(a,na,iLoc)
  
    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   
    
    INTEGER, INTENT(IN) :: iLoc
    INTEGER, INTENT(INOUT) :: na    
    INTEGER, INTENT(INOUT) :: a(:)
    
! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: ia

! ******************************************************************************
!   Start
! ******************************************************************************  

    DO ia = iLoc,na-1
      a(ia) = a(ia+1)
    END DO ! ia

    a(na) = CRAZY_VALUE_INT

    na = na - 1

! ******************************************************************************
!   End
! ******************************************************************************  
  
  END SUBROUTINE RemoveInteger









! ******************************************************************************
!
! Purpose: Delete duplicated elements of sorted integer arrays a of size na.
!
! Description: None.
!
! Input:
!   a   Array containing sorted elements, possibly containing duplicates
!   na  Number of elements in array a
!
! Output:
!   a   Array containing sorted elements without duplicates
!   nb  Number of elements in sorted array a
!
! Notes: None. 
!
! ******************************************************************************

  SUBROUTINE SimplifySortedIntegers(a,na,nb)
  
    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Parameters
! ==============================================================================   
    
    INTEGER, INTENT(IN) :: na
    INTEGER, INTENT(INOUT) :: a(na)
    INTEGER, INTENT(OUT) :: nb    
    
! ==============================================================================
!   Locals
! ==============================================================================   

    INTEGER :: ia,im
    INTEGER :: m(na)

! ******************************************************************************
!   Start
! ******************************************************************************  

    nb   = 0
    m(:) = 0
    m(1) = a(1)
    im   = 1

    DO ia = 2,na
      IF ( a(ia) /= a(ia-1) ) THEN 
        im = im + 1      
        m(im) = a(ia)
      END IF ! a
    END DO ! ia

    nb = im
    
    a(1:im) = m(1:im)    
    a(im+1:na) = 0
  
  END SUBROUTINE SimplifySortedIntegers




! ******************************************************************************
! End
! ******************************************************************************

END MODULE ModSortSearch


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModSortSearch.F90,v $
! Revision 1.18  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.15  2006/03/20 13:52:56  haselbac
! Added new routine, some clean-up
!
! Revision 1.14  2005/06/20 17:06:16  haselbac
! New routine to rm common ints from lists, put routines in proper order
!
! Revision 1.13  2005/05/17 01:12:36  haselbac
! Improved description of BinarySearchInteger, cosmetics
!
! Revision 1.12  2005/04/21 01:20:03  haselbac
! Bug fix in binary search; extended to return info if failed
!
! Revision 1.11  2004/12/29 21:03:52  haselbac
! Added new procedure, cosmetics
!
! Revision 1.10  2004/01/22 16:02:41  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.9  2003/12/04 03:28:28  haselbac
! Added routine QuickSortIntegerInteger
!
! Revision 1.8  2003/11/25 21:03:14  haselbac
! Fixed bug in MergeSortedIntegers: Was accessing m(0)
!
! Revision 1.7  2003/08/19 22:46:55  haselbac
! Added CycleList routine
!
! Revision 1.6  2003/07/22 02:02:19  haselbac
! Fixed bug in merging of sorted lists
!
! Revision 1.5  2003/03/15 17:51:38  haselbac
! Added routine to simplify sorted list of integers
!
! Revision 1.4  2003/01/28 16:47:51  haselbac
! Bug fix in MergeSortedIntegers
!
! Revision 1.3  2002/09/09 15:00:14  haselbac
! Removed some error checking, changed interface to MergeSortedIntegers
!
! Revision 1.2  2002/04/11 18:55:04  haselbac
! Added routine to sort integer array based on real key
!
! Revision 1.1  2002/03/01 17:04:24  haselbac
! Initial revision
!
! ******************************************************************************







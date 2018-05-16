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
! Purpose: Suite of routines to carry out hash table operations.
!
! Description: None.
!
! Notes: To create and use a hash table, one has to take the following steps:
!   1. Allocate a table of the appropriate size by calling CreateHashTable. 
!      This routine will automatically make the table larger in order to 
!      reduce the number of collisions.
!   2. Provide a routine which computes a key and call hashKey.
!   3. Deallocate the table by calling DestroyHashTable.
!
! ******************************************************************************
!
! $Id: RFLU_ModHashTable.F90,v 1.17 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModHashTable

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch

  IMPLICIT NONE

  PRIVATE :: RFLU_FindNearestPrime, & 
             RFLU_HashFuncPrimary, & 
             RFLU_HashFuncSecondary
             
  PUBLIC :: RFLU_CreateHashTable, & 
            RFLU_DestroyHashTable, &
            RFLU_HashBuildKey, &  
            RFLU_HashEdge, & 
            RFLU_HashFace, & 
            RFLU_HashVertex, &
            RFLU_HashVertexFancy, &
            RFLU_UnHashBFace
  
  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Public
! ==============================================================================
    
  INTEGER, PARAMETER, PUBLIC :: HASHTABLE_ENTRYSTATUS_OLD = 0, & 
                                HASHTABLE_ENTRYSTATUS_NEW = 1  
    
  INTEGER, PUBLIC :: hashTableCollisions,hashTableSize
  INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: hashTable

! ==============================================================================
! Private
! ==============================================================================
     
  INTEGER, PARAMETER, PRIVATE :: HASHTABLE_INIT = 0 ! must be <= 0      
  INTEGER, PARAMETER, PRIVATE :: NPRIMES = 48 

  INTEGER, PRIVATE :: primeNumbers(NPRIMES) = &
             (/         1,        251,        379,        509, &
                      761,       1021,       1531,       2039, &
                     3067,       4093,       6143,       8191, &
                    12289,      16381,      24571,      32749, &
                    49139,      65521,      98297,     131071, &
                   196613,     262139,     393209,     524287, &
                   786431,    1048573,    1572853,    2097143, &
                  3145721,    4194301,    6291449,    8388593, &
                 12582917,   16777213,   25165807,   33554393, &
                 50331599,   67108859,  100663261,  134217689, &
                201326549,  268435399,  402653171,  536870909, &
                805306349, 1073741789, 1610612711, 2147483647   /)
  
        
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS
  
  
  
  
! ******************************************************************************
!
! Purpose: Create hash table.
!
! Description: None.
!
! Input:
!   global      Global pointer
!   size        Size of hash table
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateHashTable(global,size)

    IMPLICIT NONE    

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: size
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,ih

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_CreateHashTable',&
  'RFLU_ModHashTable.F90')
    
! ******************************************************************************
!   Find nearest prime
! ******************************************************************************

    CALL RFLU_FindNearestPrime(2*size,hashTableSize)      

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    ALLOCATE(hashTable(hashTableSize),STAT=errorFlag)
    global%error = errorFlag       
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'hashTable') 
    END IF ! global%error        

    DO ih = 1,hashTableSize
      hashTable(ih) = HASHTABLE_INIT
    END DO ! ih
    
    hashTableCollisions = 0      

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateHashTable







! ******************************************************************************
!
! Purpose: Destroy hash table
!
! Description: None.
!
! Input:
!   global      Global pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_DestroyHashTable(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_DestroyHashTable',&
  'RFLU_ModHashTable.F90')

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************  

    DEALLOCATE(hashTable,STAT=errorFlag)
    global%error = errorFlag       
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'hashTable') 
    END IF ! global%error          

! ******************************************************************************
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyHashTable
  
  







! ******************************************************************************
!
! Purpose: Find nearest prime for hash table size.
!
! Description: None.
!
! Input:
!   size        Original size
!
! Output: 
!   primeSize   Nearest prime to original size
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_FindNearestPrime(size,primeSize)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================
      
    INTEGER, INTENT(IN) :: size
    INTEGER, INTENT(OUT) :: primeSize

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: i 

! ******************************************************************************
!   Start, find nearest prime
! ******************************************************************************

    DO i = 1,NPRIMES-1      
      IF ( size >= primeNumbers(i) .AND. size < primeNumbers(i+1) ) THEN 
        primeSize = primeNumbers(i+1)
        EXIT
      END IF ! size
    END DO ! i   

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_FindNearestPrime






! ******************************************************************************
!
! Purpose: Compute key from series of integers.
!
! Description: None.
!
! Input:
!   a           Set of integers
!   aSize       Size of set of integers
!
! Output: 
!   key         Key
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashBuildKey(a,aSize,key)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: aSize
    INTEGER, INTENT(IN) :: a(1:aSize)
    INTEGER, INTENT(OUT) :: key

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: i,term
    INTEGER, PARAMETER :: A_RAND = 31415, B_RAND = 27183

! ******************************************************************************
!   Start, compute key
! ******************************************************************************  


!    key = MOD(MOD(A_RAND*B_RAND,hashTableSize)*v1 + v2,HUGE(1)/10)
!    key = MOD(A_RAND*v1 + B_RAND*v2,hashTableSize) 

! - ABS function used to guard against negative keys from integer overflow   
!    key = ABS(ABS(A_RAND*v1) + ABS(MOD(A_RAND*B_RAND,hashTableSize)*v2))       

    term = A_RAND
    key  = 0

    DO i = 1,aSize      
      term = MOD(term*B_RAND,hashTableSize)

      IF ( term < 0 ) THEN 
        term = term + HUGE(1) + 1
      END IF ! term

      key  = MOD((term*key + a(i)),HUGE(1))

      IF ( key < 0 ) THEN 
        key = key + HUGE(1) + 1
      END IF ! key
    END DO ! i

! ******************************************************************************
!   End
! ******************************************************************************  

  END SUBROUTINE RFLU_HashBuildKey
  



! ******************************************************************************
!
! Purpose: Compute key from integer.
!
! Description: None.
!
! Input:
!   a           Integer
!
! Output: 
!   key         Key
!
! Notes: Added as separate function because passing single-element arrays lead
!   to warnings from Intel f90 compiler.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashBuildKey1(a,key)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: a
    INTEGER, INTENT(OUT) :: key

! ==============================================================================
!   Local variables
! ==============================================================================

! ******************************************************************************
!   Start, compute key
! ******************************************************************************  

    key = MOD(a,HUGE(1))

    IF ( key < 0 ) THEN 
      key = key + HUGE(1) + 1
    END IF ! key

! ******************************************************************************
!   End
! ******************************************************************************  

  END SUBROUTINE RFLU_HashBuildKey1






! ******************************************************************************
!
! Purpose: Hash edge.
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!   pGrid       Pointer to grid
!   v           Edge vertices
!
! Output:
!   edgeType    Flag indicating whether edge is new or not
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashEdge(global,key,pGrid,v,edgeType)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(IN) :: v(1:2)
    INTEGER, INTENT(OUT) :: edgeType
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: addr,collCntr,incr
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_HashEdge',&
  'RFLU_ModHashTable.F90')
    
! ******************************************************************************
!   Construct address from key based on sorted edge vertices
! ******************************************************************************
      
    CALL RFLU_HashFuncPrimary(key,addr)  
  
! ******************************************************************************
!   Insert edge into hash table
! ******************************************************************************  
  
    collCntr = 0

    DO    
      IF ( hashTable(addr) == HASHTABLE_INIT ) THEN
        IF ( pGrid%nEdgesTot == pGrid%nEdgesEst ) THEN       
          CALL ErrorStop(global,ERR_NEDGES_ESTIMATE,__LINE__)
        END IF ! pGrid%nEdgesTot    

        edgeType = EDGE_TYPE_NEW

        pGrid%nEdgesTot = pGrid%nEdgesTot + 1     
        hashTable(addr) = pGrid%nEdgesTot

        pGrid%e2v(1,pGrid%nEdgesTot) = v(1)
        pGrid%e2v(2,pGrid%nEdgesTot) = v(2)     

        EXIT
      ELSE       
        IF ( v(1) == pGrid%e2v(1,hashTable(addr)) .AND. & 
             v(2) == pGrid%e2v(2,hashTable(addr)) ) THEN
          edgeType = EDGE_TYPE_OLD        

          EXIT
        ELSE       
          hashTableCollisions = hashTableCollisions + 1

          IF ( collCntr == 0 ) THEN 
            CALL RFLU_HashFuncSecondary(key,incr)
          END IF ! collCntr

          collCntr = collCntr + 1        
          addr = MOD(addr + incr,hashTableSize)

          IF ( addr == 0 ) THEN 
            addr = 1
          END IF 
        END IF ! fv

      END IF ! hashTable
    END DO ! <empty>

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HashEdge







! ******************************************************************************
!
! Purpose: Hash face.
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!   pGrid       Pointer to grid
!   icg         Global cell number
!   ifl         Local face number
!   fv          Face vertices (only three are needed)
!   nVert       Number of vertices in face (3 = triangle, 4 = quadrilateral)
!
! Output:
!   faceType    Flag indicating whether face is new or not
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashFace(global,key,pGrid,icg,ifl,fv,nVert,faceType)
  
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: icg,ifl,key,nVert
    INTEGER, INTENT(IN) :: fv(1:3)
    INTEGER, INTENT(OUT) :: faceType
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: addr,collCntr,incr
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_HashFace',&
  'RFLU_ModHashTable.F90')
    
! ******************************************************************************
!   Construct address from key based on sorted face vertices
! ******************************************************************************
      
    CALL RFLU_HashFuncPrimary(key,addr)  

! ******************************************************************************
!   Insert face into hash table
! ******************************************************************************  

    collCntr = 0

    DO    

! ==============================================================================
!     Entry not yet occupied
! ==============================================================================  

      IF ( hashTable(addr) == HASHTABLE_INIT ) THEN     
        IF ( pGrid%nFacesTot == pGrid%nFacesEst ) THEN       
          CALL ErrorStop(global,ERR_NFACES_ESTIMATE,__LINE__)
        END IF ! pGrid%nFacesTot    

        faceType = FACE_TYPE_NEW

        pGrid%nFacesTot = pGrid%nFacesTot + 1     
        hashTable(addr) = pGrid%nFacesTot

        pGrid%f2c(1  ,pGrid%nFacesTot) = icg
        pGrid%f2c(3  ,pGrid%nFacesTot) = ifl  
        pGrid%f2c(4  ,pGrid%nFacesTot) = nVert   
        pGrid%f2v(1:3,pGrid%nFacesTot) = fv(1:3)

        EXIT

! ==============================================================================
!     Entry already occupied  
! ==============================================================================  

      ELSE 

! ------------------------------------------------------------------------------
!       Entry occupied by same face
! ------------------------------------------------------------------------------    

        IF ( fv(1) == pGrid%f2v(1,hashTable(addr)) .AND. & 
             fv(2) == pGrid%f2v(2,hashTable(addr)) .AND. & 
             fv(3) == pGrid%f2v(3,hashTable(addr)) ) THEN
          faceType = FACE_TYPE_OLD        

          pGrid%f2c(2,hashTable(addr)) = icg

          IF ( pGrid%f2c(4,hashTable(addr)) /= nVert ) THEN 
            CALL ErrorStop(global,ERR_FACE_NVERT_INVALID,__LINE__)
          END IF ! pGrid%f2c

          EXIT

! ------------------------------------------------------------------------------
!       Entry occupied by some other face: COLLISION
! ------------------------------------------------------------------------------    

        ELSE       
          hashTableCollisions = hashTableCollisions + 1

          IF ( collCntr == 0 ) THEN 
            CALL RFLU_HashFuncSecondary(key,incr)
          END IF ! collCntr

          collCntr = collCntr + 1        
          addr = MOD(addr + incr,hashTableSize)

          IF ( addr == 0 ) THEN 
            addr = 1
          END IF ! addr   
        END IF ! fv

      END IF ! hashTable
    END DO ! <empty>
  
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HashFace







! ******************************************************************************
!
! Purpose: Primary hash function: Compute address from key value.
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!
! Output:
!   addr        Address
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashFuncPrimary(key,addr)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: addr

! ******************************************************************************
!   Start, compute address
! ******************************************************************************

!      addr = 1 ! simple test, will lead to many collisions, but works
    addr = 1 + MOD(key,hashTableSize)             

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_HashFuncPrimary







! ******************************************************************************
!
! Purpose: Secondary hash function: Compute address from key value
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!
! Output:
!   addr        Address
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashFuncSecondary(key,addr)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: addr

! ******************************************************************************
!   Start, compute address
! ******************************************************************************

!      addr = 1 ! simple test, will lead to many collisions, but works
    addr = 1 + MOD(key,hashTableSize-2)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_HashFuncSecondary
    
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Hash vertex.
!
! Description: None.
!
! Input:
!   global      Pointer to global type
!   key         Key from which address is computed
!   ivg         Global boundary vertex number
!   nVert       Number of vertices
!   vert        List of vertices
!
! Output: 
!   nVert       Number of vertices
!   vert        List of vertices
!   errorFlag   Error flag (optional)
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashVertex(global,key,ivg,nVert,vert,errorFlag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: ivg,key
    INTEGER, INTENT(INOUT) :: nVert
    INTEGER, INTENT(OUT), OPTIONAL :: errorFlag
    INTEGER, DIMENSION(:), INTENT(INOUT) :: vert
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: addr,collCntr,incr
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_HashVertex',&
  'RFLU_ModHashTable.F90')

    IF ( PRESENT(errorFlag) ) THEN 
      errorFlag = ERR_NONE
    END IF ! PRESENT
    
! ******************************************************************************
!   Construct address from key
! ******************************************************************************
      
    CALL RFLU_HashFuncPrimary(key,addr)  
  
! ******************************************************************************
!   Insert vertex into hash table
! ******************************************************************************  
  
    collCntr = 0

    emptyLoop: DO    
      IF ( hashTable(addr) == HASHTABLE_INIT ) THEN
        IF ( nVert == SIZE(vert,1) ) THEN 
          IF ( .NOT. PRESENT(errorFlag) ) THEN 
            CALL ErrorStop(global,ERR_NVERT_ESTIMATE,__LINE__)
          ELSE 
            errorFlag = ERR_NONE + 1 ! Any value other than ERR_NONE
            EXIT emptyLoop
          END IF ! errorFlag
        END IF ! nVert  

        hashTable(addr) = ivg 
        nVert = nVert + 1     

        vert(nVert) = ivg

        EXIT emptyLoop
      ELSE 
        IF ( hashTable(addr) == ivg ) THEN            
          EXIT emptyLoop
        ELSE       
          hashTableCollisions = hashTableCollisions + 1

          IF ( collCntr == 0 ) THEN 
            CALL RFLU_HashFuncSecondary(key,incr)
          END IF ! collCntr

          collCntr = collCntr + 1        
          addr = MOD(addr + incr,hashTableSize)

          IF ( addr == 0 ) THEN 
            addr = 1
          END IF 
        END IF ! hashTable

      END IF ! hashTable
    END DO emptyLoop
  
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HashVertex  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Hash vertex and return some information about entry.
!
! Description: None.
!
! Input:
!   global      Pointer to global type
!   key         Key from which address is computed
!   ivg         Global boundary vertex number
!   nVert       Number of vertices
!   vert        List of vertices
!   indx        List of vertex indices
!
! Output: 
!   nVert       Number of vertices
!   vert        List of vertices
!   indx        List of vertex indices 
!   ivgStat     Status of vertex entry
!   ivgIndx     Index of vertex entry
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HashVertexFancy(global,key,ivg,nVert,vert,indx,ivgStat, & 
                                  ivgIndx)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: ivg,key
    INTEGER, INTENT(INOUT) :: nVert
    INTEGER, INTENT(OUT) :: ivgIndx,ivgStat
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indx,vert
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: addr,collCntr,incr
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_HashVertexFancy',&
  'RFLU_ModHashTable.F90')
    
! ******************************************************************************
!   Construct address from key
! ******************************************************************************
      
    CALL RFLU_HashFuncPrimary(key,addr)  
  
! ******************************************************************************
!   Insert vertex into hash table
! ******************************************************************************  
  
    collCntr = 0

    DO    
      IF ( hashTable(addr) == HASHTABLE_INIT ) THEN
        IF ( nVert == SIZE(vert,1) ) THEN 
          CALL ErrorStop(global,ERR_NVERT_ESTIMATE,__LINE__)
        END IF ! nVert  

        hashTable(addr) = ivg 
        nVert = nVert + 1     

        vert(nVert) = ivg
        
        indx(addr) = nVert

        ivgStat = HASHTABLE_ENTRYSTATUS_NEW
        ivgIndx = CRAZY_VALUE_INT

        EXIT
      ELSE 
        IF ( hashTable(addr) == ivg ) THEN
          ivgStat = HASHTABLE_ENTRYSTATUS_OLD
          ivgIndx = indx(addr)    
                      
          EXIT
        ELSE       
          hashTableCollisions = hashTableCollisions + 1

          IF ( collCntr == 0 ) THEN 
            CALL RFLU_HashFuncSecondary(key,incr)
          END IF ! collCntr

          collCntr = collCntr + 1        
          addr = MOD(addr + incr,hashTableSize)

          IF ( addr == 0 ) THEN 
            addr = 1
          END IF 
        END IF ! hashTable

      END IF ! hashTable
    END DO ! <empty>
  
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HashVertexFancy  
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Unhash boundary face to find the global cell number to which the
!   face is attached.
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!   pGrid       Pointer to grid
!   fv          Face vertices (only three are needed)
!   nVert       Number of vertices in face (3 = triangle, 4 = quadrilateral)
!   bcType      Boundary-condition type assigned to face
!
! Output:
!   icg         Global cell number to which boundary face is attached
!   ifg         Global face number
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_UnHashBFace(global,key,pGrid,fv,nVert,bcType,icg,ifg)
  
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: bcType,key,nVert
    INTEGER, INTENT(IN) :: fv(1:3)
    INTEGER, INTENT(OUT) :: icg,ifg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: addr,collCntr,fvSize,incr
 
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_UnHashBFace',&
  'RFLU_ModHashTable.F90')
    
! ******************************************************************************
!   Construct address from key based on sorted face vertices
! ******************************************************************************
   
    CALL RFLU_HashFuncPrimary(key,addr)  
  
! ******************************************************************************
!   Insert face into hash table
! ******************************************************************************  
  
    collCntr = 0

    DO  
      ifg = hashTable(addr)

! ==============================================================================
!     Entry not occupied, must have error in hash table
! ==============================================================================

      IF ( ifg == HASHTABLE_INIT ) THEN 
        CALL ErrorStop(global,ERR_HASHTABLE,__LINE__) 

! ==============================================================================
!     Entry occupied
! ==============================================================================

      ELSE       
      
! ------------------------------------------------------------------------------
!       No collision if vertices match
! ------------------------------------------------------------------------------
            
        IF ( fv(1) == pGrid%f2v(1,ifg) .AND. & 
             fv(2) == pGrid%f2v(2,ifg) .AND. & 
             fv(3) == pGrid%f2v(3,ifg) ) THEN
          icg = pGrid%f2c(1,ifg)

! ------- If second entry IS NOT exterior cell (as initialized), then have error
!         unless face is on symmetry or periodic patch, because there have 
!         virtual cells adjacent to patch

          IF ( pGrid%f2c(2,ifg) /= CELL_TYPE_EXT ) THEN
            IF ( bcType /= BC_SYMMETRY .AND. bcType /= BC_PERIODIC ) THEN  
              CALL ErrorStop(global,ERR_HASHTABLE,__LINE__) 
            END IF ! bcType
            
! ------- If second entry IS interior cell (as initialized), then only set to 
!         boundary cell type if face is not symmetry or periodic patch. This is
!         done to make sure this boundary face shows up as a VX face in the 
!         final face list which is has to to get geometry correct. In the 
!         sketch below, the face in question is marked by '+'. This face MUST 
!         show up as a VX face. If it is marked as a VB face, it does not show
!         up in the main facelist, and is then eliminated from the geometry
!         computation because any faces (actual and virtual) are ignored 
!         (this is done because the AV and VV faces associated with sype 
!         patches appear on the patch face list AND in the volume face list).
!         But by ignoring these faces, also do not get contribution of face 
!         marked below and geometry is incorrect. With the treatment below,
!         have this face appear twice also, and this treatment thus not only
!         fixes geometry computation problem, but is also more consistent with 
!         existing treatment.          
!
!                |         |         |
!                | virtual |  real   |  real
!               vx  cell  av  cell  aa  cell
!                |         |         |
!                |         |         |
!                |+++vx++++|---av----|------   <--- sype patch
!                          |         |
!                          | virtual | virtual 
!                          |  cell  vv  cell
!                          |         |
!                          |---vx----|------
!             
            
          ELSE 
            IF ( bcType /= BC_SYMMETRY .AND. bcType /= BC_PERIODIC ) THEN 
              pGrid%f2c(2,ifg) = CELL_TYPE_BND
	    ELSE 
	      IF ( icg <= pGrid%nCells ) THEN 
	        pGrid%f2c(2,ifg) = CELL_TYPE_BND
	      END IF ! icg
            END IF ! bcType
          END IF ! pGrid%f2c

          IF ( pGrid%f2c(4,ifg) /= nVert ) THEN 
            CALL ErrorStop(global,ERR_FACE_NVERT_INVALID,__LINE__)
          END IF ! pGrid%f2c        

          EXIT

! ------------------------------------------------------------------------------
!       Collision if vertices do not match, generate new key and address
! ------------------------------------------------------------------------------
            
        ELSE 
          IF ( collCntr == 0 ) THEN 
            CALL RFLU_HashFuncSecondary(key,incr)
          END IF ! collCntr

          collCntr = collCntr + 1        
          addr = MOD(addr + incr,hashTableSize)

          IF ( addr == 0 ) THEN 
            addr = 1
          END IF 
        END IF ! fv

      END IF ! hashTable
    END DO ! <empty>

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_UnHashBFace


  
  
  
  

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModHashTable


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModHashTable.F90,v $
! Revision 1.17  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.14  2006/12/15 13:24:20  haselbac
! Added new RFLU_HashBuildKey1 for single integers to avoid ifort warnings
!
! Revision 1.13  2006/04/17 19:53:37  haselbac
! Bug fix: For serial cells, opposite cell should be BND, not EXT
!
! Revision 1.12  2006/04/13 18:09:19  haselbac
! Bug fix: Treatment of virtual faces on sype patches
!
! Revision 1.11  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.10  2006/03/25 21:53:29  haselbac
! Changes because of sype patches
!
! Revision 1.9  2005/06/14 17:46:56  haselbac
! Adapted RFLU_HashVertex to return optional error flag
!
! Revision 1.8  2005/04/21 01:36:26  haselbac
! Added routine which returns info, in particular if element already hashed
!
! Revision 1.7  2004/12/04 03:31:04  haselbac
! Changed RFLU_HashBVertex so can be used for any type of vertex list
!
! Revision 1.6  2004/07/06 15:14:42  haselbac
! Added subroutines for hashing objects, cosmetics
!                                              
! Revision 1.5  2004/01/22 16:03:59  haselbac                                  
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC   and titan
!
! Revision 1.4  2003/12/04 03:28:48  haselbac                                  
! Increase size of hash table to avoid degeneracy for one case                 
!
! Revision 1.3  2002/10/08 15:49:21  haselbac                                  
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem           
!
! Revision 1.2  2002/09/09 15:08:11  haselbac                                  
! global now under regions                                                     
!
! Revision 1.1  2002/03/01 16:31:21  haselbac                                  
! Initial revision                                                             
!
! ******************************************************************************














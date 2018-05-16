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
!******************************************************************************
!
! Purpose: reads the #RVAV_COMPARISONS section of a file (until # is 
!          encountered at the end of all the comparison lists), read
!          just ONE keyword and store the associated numerical value.
!
! Description: 
!
! Input: IF_RVAV_INPUT = file number and the file RVAV_Input.inp
!
! Output: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ReadComparisonsSection.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadComparisonsSection( global ) 

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE ModParameters
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  IMPLICIT NONE

! ... parameter variables
  TYPE(t_global), POINTER :: global
                                                      
! ... loop variables
  INTEGER :: ival, n

! ... local variables
  CHARACTER(256) :: line, key

  INTEGER :: nc
  INTEGER :: nRows, nCols, errorFlag

  LOGICAL :: defined

  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: valsS1, valsS2 

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ReadComparisonsSection',&
  'RVAV_ReadComparisonsSection.F90' )

! read lines from file until # or EOF found

  defined = .false.    ! initial status
  nCols   = 12         ! this value is fixed. each stream has 12 columns

  key = "RVAV_NUMBER_OF_COMPARISONS"
  nc  = LEN_TRIM(key)

  DO
    READ(IF_RVAV_INPUT,'(A256)',err=10,end=10) line
    IF (line(1:1) == '#') EXIT

    IF (line(1:nc) == TRIM(key)) THEN

      READ(line(nc+1:256),*,err=10,end=10)  nRows
       
      globalRVAV%nComparisons = nRows
      IF (globalRVAV%nComparisons > 0) THEN
        ALLOCATE( globalRVAV%RVAVcompare(1:globalRVAV%nComparisons), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
        ALLOCATE( valsS1(nRows,nCols),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
        ALLOCATE( valsS2(nRows,nCols),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
	
        defined = .true.

        DO ival=1,globalRVAV%nComparisons
	
! ... read in the first row corresponding to stream one (S1)
          READ(IF_RVAV_INPUT,*,err=10,end=10) (valsS1(ival,n), n=1,nCols)
	  
! ... read in the second row corresponding to stream two (S2)
          READ(IF_RVAV_INPUT,*,err=10,end=10) (valsS2(ival,n), n=1,nCols)
	  
! ... assign values that are read in to RVAVcompare(ival)%<fields>
! ... NOTE: data type RVAVcompare is defined in RVAV_ModGlobal.F90

          globalRVAV%RVAVcompare(ival)%operationS1     = valsS1(ival,1)
          globalRVAV%RVAVcompare(ival)%VariableIndexS1 = valsS1(ival,2)
          globalRVAV%RVAVcompare(ival)%blockS1         = valsS1(ival,3)
          globalRVAV%RVAVcompare(ival)%ibegS1          = valsS1(ival,4)
          globalRVAV%RVAVcompare(ival)%iendS1          = valsS1(ival,5)
          globalRVAV%RVAVcompare(ival)%ijumpS1         = valsS1(ival,6)
          globalRVAV%RVAVcompare(ival)%jbegS1          = valsS1(ival,7)
          globalRVAV%RVAVcompare(ival)%jendS1          = valsS1(ival,8)
          globalRVAV%RVAVcompare(ival)%jjumpS1         = valsS1(ival,9)
          globalRVAV%RVAVcompare(ival)%kbegS1          = valsS1(ival,10)
          globalRVAV%RVAVcompare(ival)%kendS1          = valsS1(ival,11)
          globalRVAV%RVAVcompare(ival)%kjumpS1         = valsS1(ival,12)

          globalRVAV%RVAVcompare(ival)%operationS2     = valsS2(ival,1)
          globalRVAV%RVAVcompare(ival)%VariableIndexS2 = valsS2(ival,2)
          globalRVAV%RVAVcompare(ival)%blockS2         = valsS2(ival,3)
          globalRVAV%RVAVcompare(ival)%ibegS2          = valsS2(ival,4)
          globalRVAV%RVAVcompare(ival)%iendS2          = valsS2(ival,5)
          globalRVAV%RVAVcompare(ival)%ijumpS2         = valsS2(ival,6)
          globalRVAV%RVAVcompare(ival)%jbegS2          = valsS2(ival,7)
          globalRVAV%RVAVcompare(ival)%jendS2          = valsS2(ival,8)
          globalRVAV%RVAVcompare(ival)%jjumpS2         = valsS2(ival,9)
          globalRVAV%RVAVcompare(ival)%kbegS2          = valsS2(ival,10)
          globalRVAV%RVAVcompare(ival)%kendS2          = valsS2(ival,11)
          globalRVAV%RVAVcompare(ival)%kjumpS2         = valsS2(ival,12)

! ... check read sanity

          IF ( global%verbLevel/=VERBOSE_NONE ) THEN
            WRITE(STDOUT,'(/,A)')    'ReadComparisonSection'
            WRITE(STDOUT,'(A,I5)')   'nComparisons = ',globalRVAV%nComparisons
            WRITE(STDOUT,'(/,A,I5)') 'Stream1',ival
            WRITE(STDOUT,'(A,I5)')   '  operationS1 =',globalRVAV%RVAVcompare(ival)%operationS1
            WRITE(STDOUT,'(A,I5)')   '  variableIndexS1 =',globalRVAV%RVAVcompare(ival)%VariableIndexS1
            WRITE(STDOUT,'(A,I5)')   '  ibegS1  = ',globalRVAV%RVAVcompare(ival)%ibegS1
            WRITE(STDOUT,'(A,I5)')   '  iendS1  = ',globalRVAV%RVAVcompare(ival)%iendS1
            WRITE(STDOUT,'(A,I5)')   '  ijumpS1 = ',globalRVAV%RVAVcompare(ival)%ijumpS1
            WRITE(STDOUT,'(A,I5)')   '  jbegS1  = ',globalRVAV%RVAVcompare(ival)%jbegS1
            WRITE(STDOUT,'(A,I5)')   '  jendS1  = ',globalRVAV%RVAVcompare(ival)%jendS1
            WRITE(STDOUT,'(A,I5)')   '  jjumpS1 = ',globalRVAV%RVAVcompare(ival)%jjumpS1
            WRITE(STDOUT,'(A,I5)')   '  kbegS1  = ',globalRVAV%RVAVcompare(ival)%kbegS1
            WRITE(STDOUT,'(A,I5)')   '  kendS1  = ',globalRVAV%RVAVcompare(ival)%kendS1
            WRITE(STDOUT,'(A,I5)')   '  kjumpS1 = ',globalRVAV%RVAVcompare(ival)%kjumpS1
            WRITE(STDOUT,'(/,A)')    'Stream2'
            WRITE(STDOUT,'(A,I5)')   '  operationS2     =',globalRVAV%RVAVcompare(ival)%operationS2
            WRITE(STDOUT,'(A,I5)')   '  variableIndexS2 =',globalRVAV%RVAVcompare(ival)%VariableIndexS2
            WRITE(STDOUT,'(A,I5)')   '  ibegS2  = ',globalRVAV%RVAVcompare(ival)%ibegS2
            WRITE(STDOUT,'(A,I5)')   '  iendS2  = ',globalRVAV%RVAVcompare(ival)%iendS2
            WRITE(STDOUT,'(A,I5)')   '  ijumpS1 = ',globalRVAV%RVAVcompare(ival)%ijumpS2
            WRITE(STDOUT,'(A,I5)')   '  jbegS2  = ',globalRVAV%RVAVcompare(ival)%jbegS2
            WRITE(STDOUT,'(A,I5)')   '  jendS2  = ',globalRVAV%RVAVcompare(ival)%jendS2
            WRITE(STDOUT,'(A,I5)')   '  jjumpS2 = ',globalRVAV%RVAVcompare(ival)%jjumpS2
            WRITE(STDOUT,'(A,I5)')   '  kbegS2  = ',globalRVAV%RVAVcompare(ival)%kbegS2
            WRITE(STDOUT,'(A,I5)')   '  kendS2  = ',globalRVAV%RVAVcompare(ival)%kendS2
            WRITE(STDOUT,'(A,I5)')   '  kjumpS3 = ',globalRVAV%RVAVcompare(ival)%kjumpS2
            WRITE(STDOUT,'(/,A)')    'END ReadComparisonSection'
          END IF ! verbLevel
	  
        ENDDO ! ival

        DEALLOCATE( valsS1,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
        DEALLOCATE( valsS2,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
	
      ENDIF ! nComparisons
    ENDIF   ! line

  ENDDO ! <empty>

! finalize

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )

999  CONTINUE

END SUBROUTINE RVAV_ReadComparisonsSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ReadComparisonsSection.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:46:35  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.3  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************








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
! Purpose: read in user input related to the type of data in stream1, the
!          formats of the grids and solution files for the verification and 
!          validation operations
!
! Description: none.
!
! Input: user input file.
!
! Output: global% and globalRVAV%
!
! Notes: for the moment the vav subroutines are written only for rocflo
!
!******************************************************************************
!
! $Id: RVAV_ReadSectionStream2.F90,v 1.3 2008/12/06 08:45:09 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadSectionStream2( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  USE RVAV_ModParameters ! this module holds the RVAV specific parameters like
                         ! FILE_COMPUTED, FILE_EXPERIMENTAL and FILE_ANALYTICAL
  USE RVAV_ModGlobal     ! this is a RVAV specific module that holds the 
                         ! definitions of all the RVAV variables
  IMPLICIT NONE

! ... parameter variables
  TYPE(t_global), POINTER :: global
  
! ... local variables
  INTEGER            :: nVals
  INTEGER, PARAMETER :: NVALS_MAX = 5

  CHARACTER(CHRLEN) :: keys(NVALS_MAX)

  LOGICAL           :: defined(NVALS_MAX)

  REAL(RFREAL)      :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ReadSectionStream2',&
  'RVAV_ReadSectionStream2.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX
  
  keys( 1) = 'RVAV_FLOW_TYPE'
  keys( 2) = 'RVAV_FILE_TYPE'
  keys( 3) = 'RVAV_GRID_TYPE'
  keys( 4) = 'RVAV_SOLUTION_TYPE'
  keys( 5) = 'RVAV_SIMILARITY_TYPE'

  CALL ReadSection( global, IF_RVAV_INPUT,nVals,keys,vals, defined )

  IF (defined(1)) THEN
                      globalRVAV%flowTypeS2 = FLOW_STEADY
    IF (vals(1)>0.9 ) globalRVAV%flowTypeS2 = FLOW_UNSTEADY
  ENDIF

  IF (defined(2)) THEN
    IF (vals(2)>0.0  .AND. vals(2)<10.1) globalRVAV%fileTypeS2 = FILE_COMPUTED
    IF (vals(2)>10.1 .AND. vals(2)<20.1) globalRVAV%fileTypeS2 = FILE_EXPERIMENTAL
    IF (vals(2)>20.1 .AND. vals(2)<30.1) globalRVAV%fileTypeS2 = FILE_ANALYTICAL
  ENDIF

  IF (defined(3)) THEN
                                       globalRVAV%GridFormatS2 = FORMAT_ASCII
    IF (vals(3)>0.9 .AND. vals(3)<1.1) globalRVAV%GridFormatS2 = FORMAT_BINARY
    IF (vals(3) > 1.9)                 globalRVAV%GridFormatS2 = FORMAT_HDF
  ENDIF

  IF (defined(4)) THEN
                                       globalRVAV%SolutFormatS2 = FORMAT_ASCII
    IF (vals(4)>0.9 .AND. vals(4)<1.1) globalRVAV%SolutFormatS2 = FORMAT_BINARY
    IF (vals(4)>1.9)                   globalRVAV%SolutFormatS2 = FORMAT_HDF
  ENDIF

  IF ( globalRVAV%fileTypeS2       == FILE_ANALYTICAL ) THEN
    IF ( defined(5) ) THEN
                                             globalRVAV%SimilarityTypeS2 = 0
      IF (vals(5)>130.9 .AND. vals(5)<131.1) globalRVAV%SimilarityTypeS2 = RVAV_CULICK
      IF (vals(5)>131.9 .AND. vals(5)<132.1) globalRVAV%SimilarityTypeS2 = RVAV_BLASIUS
      IF (vals(5)>132.9 )                    globalRVAV%SimilarityTypeS2 = RVAV_GAMMBUMP
    END IF 
  ENDIF 
  
! .. checking if we have read the input deck correctly
  IF ( global%verbLevel/=VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(/,A)')  '  readSectionStream2'
    WRITE(STDOUT,'(A,I5)') '   FlowTypeS2    = ',globalRVAV%FlowTypeS2
    WRITE(STDOUT,'(A,I5)') '   FileTypeS2    = ',globalRVAV%FileTypeS2
    WRITE(STDOUT,'(A,I5)') '   GridFormatS2  = ',globalRVAV%GridFormatS2
    WRITE(STDOUT,'(A,I5)') '   SolutFormatS2 = ',globalRVAV%SolutFormatS2
    WRITE(STDOUT,'(A,I5)') '   SimilarityTypeS2 = ',globalRVAV%SimilarityTypeS2
    WRITE(STDOUT,'(A,I5)') '  readSectionStream2'  
  END IF ! verbLevel

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RVAV_ReadSectionStream2

!******************************************************************************
!
! RCS Revision history
!
! $Log: RVAV_ReadSectionStream2.F90,v $
! Revision 1.3  2008/12/06 08:45:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:31  fnajjar
! Initial revision after changing case
!
! Revision 1.10  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.4  2002/08/15 19:48:06  jblazek
! Implemented grid deformation capability.
!
! Revision 1.3  2002/06/19 20:27:20  f-najjar
! Included GAMM Bump Definition
!
! Revision 1.2  2002/06/15 17:46:31  f-najjar
! Include Similarity Solution String
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!
!******************************************************************************








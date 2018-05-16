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
! $Id: RVAV_ReadSectionStream1.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadSectionStream1( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  IMPLICIT NONE

! ... parameter variables
  TYPE(t_global), POINTER :: global
  
! ... local variables
  INTEGER            :: nVals
  INTEGER, PARAMETER :: NVALS_MAX = 4

  CHARACTER(CHRLEN) :: keys(NVALS_MAX)

  LOGICAL           :: defined(NVALS_MAX)

  REAL(RFREAL)      :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ReadSectiontream1',&
  'RVAV_ReadSectionStream1.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX
  
  keys(1) = 'RVAV_FLOW_TYPE'
  keys(2) = 'RVAV_FILE_TYPE'
  keys(3) = 'RVAV_GRID_TYPE'
  keys(4) = 'RVAV_SOLUTION_TYPE'

  CALL ReadSection( global, IF_RVAV_INPUT,nVals,keys,vals,defined )

  IF (defined(1)) THEN
                      globalRVAV%FlowTypeS1 = FLOW_STEADY
    IF (vals(1)>0.9 ) globalRVAV%FlowTypeS1 = FLOW_UNSTEADY
  ENDIF

  IF (defined(2)) THEN
    IF (vals(2)>0.0  .AND. vals(2)<10.1) globalRVAV%FileTypeS1 = FILE_COMPUTED
    IF (vals(2)>10.1 .AND. vals(2)<20.1) globalRVAV%FileTypeS1 = FILE_EXPERIMENTAL
    IF (vals(2)>20.1 .AND. vals(2)<30.1) globalRVAV%FileTypeS1 = FILE_ANALYTICAL
  ENDIF

  IF (defined(3)) THEN
                                       globalRVAV%GridFormatS1 = FORMAT_ASCII
    IF (vals(3)>0.9 .AND. vals(3)<1.1) globalRVAV%GridFormatS1 = FORMAT_BINARY
    IF (vals(3)>1.9)                   globalRVAV%GridFormatS1 = FORMAT_HDF
  ENDIF

  IF (defined(4)) THEN
                                       globalRVAV%SolutFormatS1 = FORMAT_ASCII
    IF (vals(4)>0.9 .AND. vals(4)<1.1) globalRVAV%SolutFormatS1 = FORMAT_BINARY
    IF (vals(4)>1.9)                   globalRVAV%SolutFormatS1 = FORMAT_HDF
  ENDIF

! .. checking if we have read the stough correctly

  IF ( global%verbLevel/=VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(/,A)')  '  readSectionStream1'
    WRITE(STDOUT,'(A,I5)') '   FlowTypeS1    = ',globalRVAV%FlowTypeS1
    WRITE(STDOUT,'(A,I5)') '   FileTypeS1    = ',globalRVAV%FileTypeS1
    WRITE(STDOUT,'(A,I5)') '   GridFormatS1  = ',globalRVAV%GridFormatS1
    WRITE(STDOUT,'(A,I5)') '   SolutFormatS1 = ',globalRVAV%SolutFormatS1
    WRITE(STDOUT,'(A,I5)') '  readSectionStream1'  
  END IF ! verbLevel

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RVAV_ReadSectionStream1

!******************************************************************************
!
! RCS Revision history
!
! $Log: RVAV_ReadSectionStream1.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:30  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.2  2002/08/15 19:48:06  jblazek
! Implemented grid deformation capability.
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!
!******************************************************************************








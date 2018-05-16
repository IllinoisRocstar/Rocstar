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
! Purpose: read File for Stream 2 Data.
!
! Description: none.
!
! Input: 
!
! Output: Memory location with data for Stream 2
!
! ISSUE: Where do you read the Analytical and Experimental Data on Nodes or Cells
! 
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ReadFileStream2.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadFileStream2 ( regionsS1, regionsS2 )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModInterfaces, ONLY : RFLO_ReadRegionTopology, RFLO_InitInputValues, &
        ReadInputFile, RFLO_DerivedInputValues, RFLO_GetDimensDummyNodes, &
        RFLO_GetDimensDummy, RFLO_GetNodeOffset, RFLO_GetCellOffset, &
        RFLO_ReadGridRegion, RFLO_ReadSolutionRegion
  USE ModMPI
  USE ModParameters
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModDataStruct
  USE RVAV_ModInterfaces, ONLY : RVAV_ReadFileStream2Analyt, RVAV_ReadFileStream2Comput, &
                                 RVAV_ReadFileStream2Experm
  IMPLICIT NONE
  
! ... parameter variables
  TYPE (t_region), POINTER :: regionsS1(:), regionsS2(:)
  
! ... local variables
  CHARACTER(CHRLEN) :: msg
  
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  global => regionsS1(1)%global
  
  CALL RegisterFunction( global, 'RVAV_ReadFileStream2',&
  'RVAV_ReadFileStream2.F90' )

! modify all global data pertinent to Rocflo for code reuse 
  
  global%flowType    = globalRVAV%flowTypeS2
  global%gridFormat  = globalRVAV%gridFormatS2
  global%solutFormat = globalRVAV%solutFormatS2
  
  IF (global%verbLevel /= VERBOSE_NONE) THEN 
    WRITE(STDOUT,'(/,A,I5)') 'FileTypeS2    = ', globalRVAV%fileTypeS2
    WRITE(STDOUT,'(A,I5)')   'FlowTypeS2    = ', globalRVAV%flowTypeS2
    WRITE(STDOUT,'(A,I5)')   'gridFormatS2  = ', globalRVAV%gridFormatS2
    WRITE(STDOUT,'(A,I5)')   'solutFormatS2 = ', globalRVAV%solutFormatS2
  ENDIF ! verbLevel

  SELECT CASE ( globalRVAV%fileTypeS2 )
    CASE (FILE_COMPUTED) 
      IF (global%verbLevel /= VERBOSE_NONE) &
        WRITE(STDOUT,'(/,A)') 'Entering Stream2 Computed Solution...'
      CALL RVAV_readFileStream2Comput( regionsS1,regionsS2 )
        
    CASE (FILE_ANALYTICAL)
      IF (global%verbLevel /= VERBOSE_NONE) &
        WRITE(STDOUT,'(/,A)') 'Entering Stream2 Analytical Solution...'
      CALL RVAV_readFileStream2Analyt( global, regionsS2 )
        
    CASE (FILE_EXPERIMENTAL)
      IF (global%verbLevel /= VERBOSE_NONE) &
        WRITE(STDOUT,'(/,A)') 'Entering Stream2 Experimental Solution...'
      CALL RVAV_readFileStream2Experm( global, regionsS2 )

    CASE DEFAULT
      WRITE(STDOUT,'(/,A,/)') 'Current Data Structure for Stream 2 Only supports Computed, Analytical or Experimental Results.'
      STOP          
  END SELECT ! fileTypeS2                   

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global ) 

END SUBROUTINE RVAV_readFileStream2

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ReadFileStream2.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:28  fnajjar
! Initial revision after changing case
!
! Revision 1.11  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.6  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.5  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.4  2002/07/16 22:32:51  f-najjar
! Cleanup of ReadFileStream2
!
! Revision 1.3  2002/06/24 15:51:37  f-najjar
! Added IO for high Verbosity
!
! Revision 1.2  2002/06/15 17:43:44  f-najjar
! Grid & Solution for all regions in one file
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************








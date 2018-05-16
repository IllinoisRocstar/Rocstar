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
! Purpose: write dimension file for Rocflu run.
!
! Description: none.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: TFLU_WriteFluDimens.F90,v 1.4 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE WriteFluDimens( global )

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModMPI
    
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: igPatch
   
! ... local variables
  CHARACTER(CHRLEN) :: iFileName, sectionString

  INTEGER :: errorFlag, iFile
  REAL(RFREAL) :: ratioMax2Tot
  
! ******************************************************************************

  CALL RegisterFunction(global,'WriteFluDimens',&
  'TFLU_WriteFluDimens.F90')

! start ------------------------------------------------------------------------

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Rocflu dimension file...'
  END IF 

  iFile = IF_DIMS

  WRITE(iFileName,'(A,I5.5)') &
        TRIM(global%outDir)//TRIM(global%casename)//'.dim_',0

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  global%error = errorFlag        

  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error

! header and general information -----------------------------------------------

  ratioMax2Tot = 1.2_RFREAL ! Hard-coded for now

  IF ( global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel 

  sectionString = '# ROCFLU dimensions file'  
  WRITE(iFile,'(A)') TRIM(sectionString)  

  sectionString = '# Vertices'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(3(I8))') global%tofluNVerts, global%tofluNVerts, &
                         INT(ratioMax2Tot*global%tofluNVerts)

  sectionString = '# Cells'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(3(I8))') global%tofluNHexs, global%tofluNHexs, &
                         INT(ratioMax2Tot*global%tofluNHexs)

  sectionString = '# Tetrahedra'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(3(I8))') 0, 0, 0

  sectionString = '# Hexahedra'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(3(I8))') global%tofluNHexs, global%tofluNHexs, &
                         INT(ratioMax2Tot*global%tofluNHexs)

  sectionString = '# Prisms'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(3(I8))') 0, 0, 0

  sectionString = '# Pyramids'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(3(I8))') 0, 0, 0

!  sectionString = '# Faces'
!  WRITE(iFile,'(A)') TRIM(sectionString) 
!  WRITE(iFile,'(2(I8))') global%tofluNFaces, global%tofluNFaces

!  sectionString = '# Edges'
!  WRITE(iFile,'(A)') TRIM(sectionString) 
!  WRITE(iFile,'(2(I8))') global%tofluNEdges, global%tofluNEdges

  sectionString = '# Patches'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(2(I8))') global%tofluNPatches, global%tofluNPatches

  DO igPatch = 1,global%tofluNPatches

    WRITE(iFile,'(5(I8))') igPatch, 0, 0, & 
                           global%tofluNbFaces(igPatch), & 
                           global%tofluNbFaces(igPatch)
  END DO ! iPatch

  sectionString = '# Borders'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(1(I8))') 0

! end marker -------------------------------------------------------------------

  IF ( global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF

  sectionString = '# End'
  WRITE(iFile,'(A)') TRIM(sectionString) 

! close file -------------------------------------------------------------------

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF

  IF ( global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Rocflu dimension file done.'
  END IF

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )
  
END SUBROUTINE WriteFluDimens

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_WriteFluDimens.F90,v $
! Revision 1.4  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/12/21 19:18:18  wasistho
! modified to adapt changes in Rocflu
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:58:20  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/18 02:15:12  wasistho
! added new routines to create dimension file
!
!
! ******************************************************************************








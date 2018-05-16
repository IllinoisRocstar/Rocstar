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
! Purpose: write unstructured grid file in Rocflu format.
!
! Description: none.
!
! Notes: currently only ASCII Rocflu format is supported.
!
! ******************************************************************************
!
! $Id: TFLU_WriteFluGrid.F90,v 1.4 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE WriteFluGrid( global )

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModMPI
    
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: iPatch, i, j, k
   
! ... local variables
  CHARACTER(CHRLEN) :: iFileName, sectionString

  INTEGER :: errorFlag, iFile
  INTEGER :: nTets, nPris, nPyrs, vertFlag, hexFlag
  INTEGER :: nBTris, nBQuads, nBVerts, bvFlag
  
! ******************************************************************************

  CALL RegisterFunction(global,'WriteFluGrid',&
  'TFLU_WriteFluGrid.F90')

! start ------------------------------------------------------------------------

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII Rocflu grid file...'
  END IF 

  iFile = IF_GRID

  WRITE(iFileName,'(A,I5.5)') &
        TRIM(global%outDir)//TRIM(global%casename)//'.grda_',0

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  global%error = errorFlag        

  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error

! header and general information -----------------------------------------------

  IF ( global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel 

  sectionString = '# ROCFLU grid file'  
  WRITE(iFile,'(A)') TRIM(sectionString)  

  sectionString = '# Precision and range'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

  sectionString = '# Physical time'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(E23.16)') global%currentTime 

! dimensions -------------------------------------------------------------------

  nTets = 0
  nPris = 0
  nPyrs = 0
  sectionString = '# Dimensions'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(5(I8))') global%tofluNVerts,nTets,global%tofluNHexs, &
                         nPris,nPyrs

! coordinates ------------------------------------------------------------------

  IF ( global%verbLevel > VERBOSE_NONE ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
  END IF

  sectionString = '# Coordinates'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  DO i = 1,3
    WRITE(iFile,'(5(E23.16))') (global%tofluXyz(i,j),j=1,global%tofluNVerts)
  END DO

!  vertFlag = 1
!  WRITE(iFile,'(10(I8))') (vertFlag,j=1,global%tofluNVerts)

! connectivity -----------------------------------------------------------------

  IF ( global%tofluNHexs > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...'
    END IF

    sectionString = '# Hexahedra'
    WRITE(iFile,'(A)') TRIM(sectionString)  
    DO i = 1,8
      WRITE(iFile,'(10(I8))') (global%tofluHex2v(i,j),j=1,global%tofluNHexs)
    END DO ! i

!    hexFlag = 1
!    WRITE(iFile,'(10(I8))') (hexFlag,j=1,global%tofluNHexs)    
  END IF ! nHexs

! boundary information ---------------------------------------------------------

  IF ( global%verbLevel > VERBOSE_LOW ) THEN   
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
  END IF

  sectionString = '# Boundaries'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(I8)') global%tofluNPatches

  DO iPatch = 1, global%tofluNPatches
    nBTris  = 0
    nBQuads = global%tofluNbFaces(iPatch)
!    nBVerts = global%tofluNbVerts(iPatch)

!    WRITE(iFile,'(3(I8))') nBTris, nBQuads, nBVerts
    WRITE(iFile,'(3(I8))') nBTris, nBQuads

    IF ( nBQuads > 0 ) THEN 
      DO j = 1,4
        WRITE(iFile,'(10(I8))') (global%tofluQuad2v(j,k,iPatch),k=1,nBQuads)
      END DO ! j      
    END IF ! bound 

!    WRITE(iFile,'(10(I8))') (global%tofluBLoc2g(k,iPatch),k=1,nBVerts)  

!    bvFlag = 1
!    WRITE(iFile,'(10(I8))') (bvFlag,k=1,nBVerts)           
  ENDDO ! iPatch

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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII Rocflu grid file done.'
  END IF

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )
  
END SUBROUTINE WriteFluGrid

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_WriteFluGrid.F90,v $
! Revision 1.4  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/12/21 19:18:30  wasistho
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
! Revision 1.1.1.1  2004/08/17 01:41:39  wasistho
! initial checkin
!
!
! ******************************************************************************








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
! Purpose: write cell-mapping file for Rocflu run.
!
! Description: none.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE WriteFluCellMap( global )

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModMPI
    
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: j
   
! ... local variables
  CHARACTER(CHRLEN) :: iFileName, sectionString

  INTEGER :: errorFlag, iFile, nTets, nPris, nPyrs
  
! ******************************************************************************

  CALL RegisterFunction(global,'WriteFluCellMap',&
  'TFLU_WriteFluCellMap.F90')

! start ------------------------------------------------------------------------

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Rocflu cell-mapping file...'
  END IF 

  iFile = IF_DIMS

  WRITE(iFileName,'(A,I5.5)') &
        TRIM(global%outDir)//TRIM(global%casename)//'.cmp_',0

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  global%error = errorFlag        

  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error

! header and general information -----------------------------------------------

  IF ( global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel 

  sectionString = '# ROCFLU cell mapping file'  
  WRITE(iFile,'(A)') TRIM(sectionString)  

  nTets = 0
  nPris = 0
  nPyrs = 0
  sectionString = '# Dimensions'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(4(I8))') nTets,global%tofluNHexs,nPris,nPyrs

  sectionString = '# Hexahedra'
  WRITE(iFile,'(A)') TRIM(sectionString) 
  WRITE(iFile,'(10(I8))') (j,j=1,global%tofluNHexs)

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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Rocflu cell-mapping file done.'
  END IF

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )
  
END SUBROUTINE WriteFluCellMap

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_WriteFluCellMap.F90,v $
! Revision 1.3  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/12/21 22:39:03  wasistho
! initial import
!
!
! ******************************************************************************








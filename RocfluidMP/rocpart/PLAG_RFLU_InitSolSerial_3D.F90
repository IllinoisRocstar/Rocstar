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
! Purpose: Initialize particle solution in serial region in 3d.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolSerial_3D.F90,v 1.3 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2005-2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolSerial_3D(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag
  USE ModParameters

  USE RFLU_ModFaceList, ONLY: RFLU_CreateCell2FaceList, & 
                              RFLU_BuildCell2FaceList, & 
                              RFLU_DestroyCell2FaceList    
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell
  USE RFLU_ModOctree

  USE PLAG_ModParameters    
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ccSize,errorFlag,icg,icl,iPcl,iVar
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  REAL(RFREAL) :: delFrac,xDel,xMax,xMin,xPcl,yDel,yMax,yMin,yPcl,zDel,zMax, &
                  zMin,zPcl
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolSerial_3D.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolSerial_3D', &
                        'PLAG_RFLU_InitSolSerial_3D.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing particle solution '// &
                             'for serial region in 3d...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and values
! ******************************************************************************
  
  pGrid => pRegion%grid
  
  pPlag => pRegion%plag  
    
  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

  delFrac = 0.01_RFREAL

! ******************************************************************************
! Get (slightly enlarged) grid bounding box
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

  xDel = xMax - xMin
  yDel = yMax - yMin
  zDel = zMax - zMin    

  xMin = xMin - delFrac*xDel
  xMax = xMax + delFrac*xDel
  yMin = yMin - delFrac*yDel
  yMax = yMax + delFrac*yDel
  zMin = zMin - delFrac*zDel
  zMax = zMax + delFrac*zDel    
       
! ******************************************************************************
! Build Octree and cell-to-face list
! ******************************************************************************

  CALL RFLU_CreateOctree(global,pGrid%nCells)

  CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(YCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(ZCOORD,1:pGrid%nCells), & 
                        xMin,xMax,yMin,yMax,zMin,zMax)

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error  

  CALL RFLU_CreateCell2FaceList(pRegion)
  CALL RFLU_BuildCell2FaceList(pRegion)

! ******************************************************************************
! Loop over particles in serial region
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls    

! ==============================================================================
!   Get particle location and get nearest cells from Octree
! ==============================================================================

    xPcl = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yPcl = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zPcl = pPlag%cv(CV_PLAG_ZPOS,iPcl)

    CALL RFLU_QueryOctree(xPcl,yPcl,zPcl,ccSize,cc)

! ==============================================================================
!   Perform in-cell test on cells
! ==============================================================================

    cellLoop: DO icl = 1,ccSize   
      icg = cc(icl)

      IF ( RFLU_ICT_TestInCell(pRegion,xPcl,yPcl,zPcl,icg) .EQV. .TRUE. ) THEN
	pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
	pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal                         

        EXIT cellLoop
      END IF ! RFLU_ICT_TestInCell
    END DO cellLoop
  END DO ! iPclSerial
    
! ******************************************************************************
! Destroy Octree and cell-to-face list
! ******************************************************************************

  CALL RFLU_DestroyOctree(global)  

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error   

  CALL RFLU_DestroyCell2FaceList(pRegion)        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution for serial region in 3d done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolSerial_3D

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolSerial_3D.F90,v $
! Revision 1.3  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/03/15 21:58:29  haselbac
! Initial revision
!
! Revision 1.4  2007/03/09 20:32:37  fnajjar
! Fixed incorrect registration name
!
! Revision 1.3  2006/05/22 15:33:09  fnajjar
! Fixed bug for uninitialized delFrac
!
! Revision 1.2  2006/05/05 18:41:10  haselbac
! Fix from last ci
!
! Revision 1.1  2006/05/05 17:27:57  haselbac
! Initial revision
!
! ******************************************************************************








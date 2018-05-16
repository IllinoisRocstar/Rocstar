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
! Purpose: Initialize particle solution in a region by copying data from 
!   serial region.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to parallel region
!   pRegionSerial       Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolFromSerial.F90,v 1.7 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2005-2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolFromSerial(pRegion,pRegionSerial)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag
  USE ModParameters
  USE ModSortSearch

  USE RFLU_ModRenumberings

  USE PLAG_ModParameters    
  
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: sortFlag
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,icgs,icgsMax,icgsMin,icgsNzPclLow,icgsNzPclMax, &
             icgsNzPclMin,icgsNzPclUpp,icl,iLoc,iPclPerCellCSR, &
             iPclPerCellCSRLow,iPclPerCellCSRUpp,iPclSerial,iSc2pc,iSc2pcLow, &
             iSc2pcUpp,iVar,j
  REAL(RFREAL) :: delFrac,xDel,xMax,xMin,xPcl,xPclMax,xPclMin,yDel,yMax,yMin, &
                  yPcl,yPclMax,yPclMin,zDel,zMax,zMin,zPcl,zPclMax,zPclMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag,pPlagSerial

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolFromSerial.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolFromSerial', &
                        'PLAG_RFLU_InitSolFromSerial.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing particle solution '// &
                             'from serial region...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and values
! ******************************************************************************
  
  pGrid       => pRegion%grid
  pPlag       => pRegion%plag  
  pPlagSerial => pRegionSerial%plag

  sortFlag = .FALSE.
    
  pPlag%nPcls = 0

  delFrac = 0.05_RFREAL
 
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
! Get particle bounding box 
! ******************************************************************************
  
  xPclMin =  HUGE(1.0_RFREAL)
  xPclMax = -HUGE(1.0_RFREAL)
  yPclMin =  HUGE(1.0_RFREAL)
  yPclMax = -HUGE(1.0_RFREAL)
  zPclMin =  HUGE(1.0_RFREAL)
  zPclMax = -HUGE(1.0_RFREAL)    
  
  DO iPclSerial = 1,pPlagSerial%nPcls
    xPclMin = MIN(xPclMin,pPlagSerial%cv(CV_PLAG_XPOS,iPclSerial)) 
    xPclMax = MAX(xPclMax,pPlagSerial%cv(CV_PLAG_XPOS,iPclSerial))
    yPclMin = MIN(yPclMin,pPlagSerial%cv(CV_PLAG_YPOS,iPclSerial)) 
    yPclMax = MAX(yPclMax,pPlagSerial%cv(CV_PLAG_YPOS,iPclSerial))
    zPclMin = MIN(zPclMin,pPlagSerial%cv(CV_PLAG_ZPOS,iPclSerial)) 
    zPclMax = MAX(zPclMax,pPlagSerial%cv(CV_PLAG_ZPOS,iPclSerial))            
  END DO ! iPclSerial  

! ******************************************************************************  
! If particle bounding box contained at least partially in grid bounding box,
! search for particles
! ******************************************************************************

  IF ( ((xPclMin < xMax) .AND. (yPclMin < yMax) .AND. (zPclMin < zMax)) .AND. & 
       ((xPclMax > xMin) .AND. (yPclMax > yMin) .AND. (zPclMax > zMin)) ) THEN
       
! ==============================================================================
!   Read Pxxx2Sxxx maps so can build sc2pc map which will be used to copy 
!   particle data, find range of serial indices associated with cells in this
!   region, and range of serial cells with non-zero particles
! ==============================================================================
  
    CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)
    CALL RFLU_RNMB_CreatePC2SCMap(pRegion)
    CALL RFLU_RNMB_CreatePV2SVMap(pRegion)

    CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

    CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)
    CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)

    CALL RFLU_RNMB_BuildSC2PCMap(pRegion,sortFlag)

    icgsNzPclMin = pPlagSerial%icgNzPcl(1)
    icgsNzPclMax = pPlagSerial%icgNzPcl(pPlagSerial%nCellsNzPcl)

    icgsMin = pGrid%sc2pc(1,1)
    icgsMax = pGrid%sc2pc(1,pGrid%nCells)

! ==============================================================================
!   If cell ranges overlap, then will have particles in this region. Then find
!   the lower and upper indices in the sc2pc and icgNzPcl maps to reduce 
!   looping and searching.
! ==============================================================================
  
    IF ( (icgsNzPclMin < icgsMax) .AND. (icgsNzPclMax > icgsMin) ) THEN 
      CALL BinarySearchInteger(pGrid%sc2pc(1,1:pGrid%nCells),pGrid%nCells, &
                               icgsNzPclMin,iLoc,j)
      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
        iSc2pcLow = MAX(1,MIN(j,pGrid%nCells))
      ELSE 
        iSc2pcLow = iLoc
      END IF ! iLoc

      CALL BinarySearchInteger(pGrid%sc2pc(1,1:pGrid%nCells),pGrid%nCells, &
                               icgsNzPclMax,iLoc,j)
      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
        iSc2pcUpp = MAX(1,MIN(j,pGrid%nCells))
      ELSE 
        iSc2pcUpp = iLoc
      END IF ! iLoc

      CALL BinarySearchInteger(pPlagSerial%icgNzPcl(1:pPlagSerial%nCellsNzPcl), & 
                               pPlagSerial%nCellsNzPcl,icgsMin,iLoc,j)
      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
        icgsNzPclLow = MAX(1,MIN(j,pGrid%nCells))
      ELSE 
        icgsNzPclLow = iLoc
      END IF ! iLoc

      CALL BinarySearchInteger(pPlagSerial%icgNzPcl(1:pPlagSerial%nCellsNzPcl), & 
                               pPlagSerial%nCellsNzPcl,icgsMax,iLoc,j)
      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
        icgsNzPclUpp = MAX(1,MIN(j,pGrid%nCells))
      ELSE 
        icgsNzPclUpp = iLoc
      END IF ! iLoc

! ------------------------------------------------------------------------------
!     Loop over range of cells in this region which lie in range of cells in 
!     serial region with non-zero particles. If cell actually has non-zero
!     particles, then copy all of them to this region. NOTE check whether 
!     region index of serial particle is not crazy value anymore, which 
!     indicates that serial particle was already assigned to another region, 
!     and hence indicates an erroneous double assignment.
! ------------------------------------------------------------------------------

      DO iSc2Pc = iSc2pcLow,iSc2pcUpp
        icgs = pGrid%sc2pc(1,iSc2pc)

        CALL BinarySearchInteger(pPlagSerial%icgNzPcl(icgsNzPclLow:icgsNzPclUpp), &
                                 icgsNzPclUpp-icgsNzPclLow+1,icgs,iLoc)

        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN  
          iLoc = iLoc + icgsNzPclLow - 1 

          IF ( iLoc > 1 ) THEN 
            iPclPerCellCSRLow = pPlagSerial%iPclPerCellCSRInfo(iLoc-1)+1
          ELSE 
            iPclPerCellCSRLow = 1 
          END IF ! iLoc

          iPclPerCellCSRUpp = pPlagSerial%iPclPerCellCSRInfo(iLoc)

          DO iPclPerCellCSR = iPclPerCellCSRLow,iPclPerCellCSRUpp
            iPclSerial = pPlagSerial%iPclPerCellCSR(iPclPerCellCSR)

            pPlag%nPcls = pPlag%nPcls + 1

            DO iVar = 1,pRegion%plag%nCv
              pPlag%cv(iVar,pPlag%nPcls) = pPlagSerial%cv(iVar,iPclSerial)
            END DO ! iVar

            DO iVar = 1,pRegion%plag%nArv
              pPlag%arv(iVar,pPlag%nPcls) = pPlagSerial%arv(iVar,iPclSerial)
            END DO ! iVar

            DO iVar = 1,pRegion%plag%nAiv
              pPlag%aiv(iVar,pPlag%nPcls) = pPlagSerial%aiv(iVar,iPclSerial)
            END DO ! iVar
  
            pPlag%aiv(AIV_PLAG_ICELLS,pPlag%nPcls) = pGrid%sc2pc(2,iSc2pc)
            pPlag%aiv(AIV_PLAG_REGINI,pPlag%nPcls) = pRegion%iRegionGlobal 
  
            IF ( pPlagSerial%aiv(AIV_PLAG_REGINI,iPclSerial) /= 0 ) THEN 
              CALL ErrorStop(global,ERR_PLAG_DSTR_INVALID,__LINE__)
            ELSE 
              pPlagSerial%aiv(AIV_PLAG_REGINI,iPclSerial) = & 
                pRegion%iRegionGlobal
            END IF ! pPlagSerial%aiv
          END DO ! iPclSerial
        END IF ! iloc
      END DO ! iSc2Pc 
    END IF ! icgsNzMin

    CALL RFLU_RNMB_DestroySC2PCMap(pRegion)
  END IF ! xPcl

! ******************************************************************************
! Write info
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Number of particles:', &
                                   pPlag%nPcls
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing particle solution '// &
                             'from serial region done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolFromSerial

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolFromSerial.F90,v $
! Revision 1.7  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/03/27 00:21:33  haselbac
! Adaptation to new initialization
!
! Revision 1.4  2006/05/22 15:33:09  fnajjar
! Fixed bug for uninitialized delFrac
!
! Revision 1.3  2006/05/05 17:36:24  haselbac
! Changed so do not need access to serial grid
!
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2005/05/18 22:27:45  fnajjar
! Initial revision
!
! ******************************************************************************








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
! Purpose: Search algorithm for cell indices.
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%plag%aiv
!         region%levels%plag%aivOld
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_GetCellIndices.F90,v 1.4 2009/10/26 00:19:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_getCellIndices( region, iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE PLAG_ModInterfaces, ONLY : PLAG_inCellTest, PLAG_InCellTestRobust
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iReg

! ... loop variables
  INTEGER :: i, iCellLev, iLevels, iPcls, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iCOff, ijCOff, ijkC, iNOff, ijNOff, iLev,                        &
             ijkNR, ijkNRI, ijkNRJ,ijkNRK, iHigh, iLow, iSkip, ipcbeg,ipcend, &
             jpcbeg,jpcend, jHigh, jLow, jSkip,                               &
             kpcbeg,kpcend, kHigh, kLow, kSkip,                               &   
             nCont, nCellLev, nCellLevMax, nLevels, nPcls 
  INTEGER,          DIMENSION(4)   :: indexCurr, indexNew, indexSearch
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld
  
  LOGICAL :: cellLocate, cellLocateRobust
  
  REAL(RFREAL) :: massL,tauLR,diamL

  REAL(RFREAL),          DIMENSION(3)   :: posPlag
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pCvOld
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pDv
  
  TYPE(t_plag),   POINTER :: pPlag  
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion  
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_GetCellIndices.F90,v $ $Revision: 1.4 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_getCellIndices',&
  'PLAG_GetCellIndices.F90' )

! Get dimensions --------------------------------------------------------------

  iLev  = region%currLevel

  nPcls = region%levels(iLev)%plag%nPcls 
  nCellLevMax = 1
  nLevels     = 8
  
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  
! Set pointers ----------------------------------------------------------------
   
  pPlag   => region%levels(iLev)%plag

  pAiv    => pPlag%aiv 
  pAivOld => pPlag%aivOld
  pCv     => pPlag%cv
  pCvOld  => pPlag%cvOld
  pDv     => pPlag%dv
 
! Loop over Lagrangian particles ----------------------------------------------

  DO iPcls = 1, nPcls      

! - Load current index --------------------------------------------------------
    
    indexCurr(1) = pAivOld(AIV_PLAG_ICELLS,iPcls) 
    indexCurr(2) = pAivOld(AIV_PLAG_INDEXI,iPcls) 
    indexCurr(3) = pAivOld(AIV_PLAG_INDEXJ,iPcls) 
    indexCurr(4) = pAivOld(AIV_PLAG_INDEXK,iPcls) 

! - Load particle positions ---------------------------------------------------

    posPlag(XCOORD) = pCv(CV_PLAG_XPOS,iPcls)
    posPlag(YCOORD) = pCv(CV_PLAG_YPOS,iPcls)
    posPlag(ZCOORD) = pCv(CV_PLAG_ZPOS,iPcls)
 
! - Set particle status -------------------------------------------------------

       pAiv(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_KEEP
    pAivOld(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_KEEP

! - Search current and surrounding cells --------------------------------------
    
    nCellLev = 1
    
    DO WHILE(nCellLev <= nCellLevMax)
    
      DO iLevels = 1, nLevels 
      
        SELECT CASE(iLevels)
        
! -- current cell -------------------------------------------------------------
          CASE(1)
            iLow  = indexCurr(2)
            iHigh = indexCurr(2)
            iSkip = 1
            
            jLow  = indexCurr(3)
            jHigh = indexCurr(3)
            jSkip = 1 
            
            kLow  = indexCurr(4)
            kHigh = indexCurr(4)
            kSkip = 1

! -- i-cell shift -------------------------------------------------------------
          CASE(2)
            iLow  = indexCurr(2)-nCellLev
            iHigh = indexCurr(2)+nCellLev
            iSkip = nCellLev+1
            
            jLow  = indexCurr(3)
            jHigh = indexCurr(3)
            jSkip = 1 
            
            kLow  = indexCurr(4)
            kHigh = indexCurr(4)
            kSkip = 1

! -- j-cell shift -------------------------------------------------------------
          CASE(3)
            iLow  = indexCurr(2)
            iHigh = indexCurr(2)
            iSkip = 1
            
            jLow  = indexCurr(3)-nCellLev
            jHigh = indexCurr(3)+nCellLev
            jSkip = nCellLev+1
            
            kLow  = indexCurr(4)
            kHigh = indexCurr(4)
            kSkip = 1
 
! -- k-cell shift -------------------------------------------------------------
          CASE(4)
            iLow  = indexCurr(2)
            iHigh = indexCurr(2)
            iSkip = 1
            
            jLow  = indexCurr(3)
            jHigh = indexCurr(3)
            jSkip = 1
            
            kLow  = indexCurr(4)-nCellLev
            kHigh = indexCurr(4)+nCellLev
            kSkip = nCellLev+1
                           
! -- ij edge cell shift -------------------------------------------------------------
          CASE(5)
            iLow  = indexCurr(2)-nCellLev
            iHigh = indexCurr(2)+nCellLev
            iSkip = nCellLev+1
    
            jLow  = indexCurr(3)-nCellLev
            jHigh = indexCurr(3)+nCellLev
            jSkip = nCellLev+1
     
            kLow  = indexCurr(4)
            kHigh = indexCurr(4)
            kSkip = 1
    
! -- ik edge cell shift -------------------------------------------------------
          CASE(6)
            iLow  = indexCurr(2)-nCellLev
            iHigh = indexCurr(2)+nCellLev
            iSkip = nCellLev+1
            
            jLow  = indexCurr(3)
            jHigh = indexCurr(3)
            jSkip = 1
            
            kLow  = indexCurr(4)-nCellLev
            kHigh = indexCurr(4)+nCellLev
            kSkip = nCellLev+1

! -- jk edge cell shift -------------------------------------------------------
          CASE(7) 
            iLow  = indexCurr(2)
            iHigh = indexCurr(2)
            iSkip = 1
                           
            jLow  = indexCurr(3)-nCellLev
            jHigh = indexCurr(3)+nCellLev
            jSkip = nCellLev+1
     
            kLow  = indexCurr(4)-nCellLev
            kHigh = indexCurr(4)+nCellLev
            kSkip = nCellLev+1

! -- ijk corner cell shift ----------------------------------------------------  
          CASE(8) 
            iLow  = indexCurr(2)-nCellLev
            iHigh = indexCurr(2)+nCellLev
            iSkip = nCellLev+1
     
            jLow  = indexCurr(3)-nCellLev
            jHigh = indexCurr(3)+nCellLev
            jSkip = nCellLev+1
    
            kLow  = indexCurr(4)-nCellLev
            kHigh = indexCurr(4)+nCellLev
            kSkip = nCellLev+1
            
        END SELECT !iLevels 

! -- Perform search -----------------------------------------------------------
           
        DO k=kLow,kHigh,kSkip
        DO j=jLow,jHigh,jSkip
        DO i=iLow,iHigh,iSkip
          ijkC   = IndIJK(i,  j  ,k  ,iCOff,ijCOff)
          ijkNR  = IndIJK(i,  j  ,k  ,iNOff,ijNOff)              
          ijkNRI = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkNRJ = IndIJK(i,j+1  ,k  ,iNOff,ijNOff)
          ijkNRK = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)

          indexSearch(1) = ijkC
          indexSearch(2) = i
          indexSearch(3) = j
          indexSearch(4) = k
  
          CALL PLAG_inCellTest(region, posPlag, indexSearch, &
                               ijkNR,ijkNRI,ijkNRJ,ijkNRK,   &
                               indexNew,cellLocate)

          IF ( cellLocate .EQV. .TRUE. ) GOTO 999
        
        END DO ! i 
        END DO ! j
        END DO ! k

      END DO ! iLevels  
        
      nCellLev = nCellLev+1

    END DO ! nCellLev

! - Trap error if unable to locate --------------------------------------------
     
    IF ( cellLocate .EQV. .FALSE. ) THEN       
      massL = SUM( pPlag%cv(pPlag%cvPlagMass(:),iPcls) )
      diamL = pDv(DV_PLAG_DIAM,iPcls)

      IF ( diamL > 1.0E-14_RFREAL ) THEN
        tauLR = 3.0_RFREAL*global%pi*pPlag%tv(TV_PLAG_MUELMIXT,iPcls)*diamL/massL
      ELSE
        tauLR = 0.0_RFREAL
      ENDIF

      WRITE(STDOUT,'(A)') &
       '##### Rocpart Warning: PLAG_inCellTest Unable to Locate Cell #####'   
      WRITE(STDOUT,1010) 'Default Search Failed to Locate Cell for Particle',&
       global%currentTime, iReg, iPcls,&
       pAivOld(AIV_PLAG_PIDINI,iPcls),          &
       pAivOld(AIV_PLAG_REGINI,iPcls),          &
       pAivOld(AIV_PLAG_ICELLS,iPcls),          &
       pAivOld(AIV_PLAG_INDEXI,iPcls),          &
       pAivOld(AIV_PLAG_INDEXJ,iPcls),          &
       pAivOld(AIV_PLAG_INDEXK,iPcls),          &
       pCvOld(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls), &
       posPlag(1:3),                            &
       pDv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls),    &
       pDv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls),&
       pDv(DV_PLAG_DIAM,iPcls),massL,tauLR
    
! -- Invoke robust cell search for particles potentially located --------------
!      in dummy cells ---------------------------------------------------------

      CALL PLAG_inCellTestRobust( region, posPlag,indexCurr, &
                                  indexNew,cellLocateRobust  )

! TEMPORARY
! Delete Particle from data structure failing robust cell search
!  due to clobbered geometry datastructure. Relevant for complex geometries.
!
      WRITE(STDOUT,'(A)') &
       '##### Rocpart Warning: Deleting Particle Following Robust Cell Search #####'
         pAiv(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_DELETE
      pAivOld(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_DELETE
! END TEMPORARY

      IF ( cellLocateRobust .EQV. .TRUE. ) GOTO 999
      
        WRITE(STDOUT,'(A)') &
         '##### Rocpart Error: PLAG_inCellTestRobust Unable to Locate Cell with Robust Search #####'   
        WRITE(STDOUT,1010) 'Unable to Locate Cell for Particle',&
        global%currentTime, iReg, iPcls,&
        pAivOld(AIV_PLAG_PIDINI,iPcls),          &
        pAivOld(AIV_PLAG_REGINI,iPcls),          &
        pAivOld(AIV_PLAG_ICELLS,iPcls),          &
        pAivOld(AIV_PLAG_INDEXI,iPcls),          &
        pAivOld(AIV_PLAG_INDEXJ,iPcls),          &
        pAivOld(AIV_PLAG_INDEXK,iPcls),          &
        pCvOld(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls), &
        posPlag(1:3),                            &
        pDv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls),    &
        pDv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls)
      
!        CALL ErrorStop( global,ERR_PLAG_CELLINDEX,__LINE__ )

! TEMPORARY
! Delete Particle from data structure for testing
         pAiv(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_DELETE 
      pAivOld(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_DELETE 
! END TEMPORARY

    END IF ! cellLocate 

! - Load new cell indices -----------------------------------------------------
         
999 CONTINUE

!  Dont load aiv into aivOld till end of RK-stage since 
!   search algorithm is modulated by RK time-stepping that
!   is based on cvOld
!
!    pAivOld(AIV_PLAG_ICELLS,iPcls) = pAiv(AIV_PLAG_ICELLS,iPcls)
!    pAivOld(AIV_PLAG_INDEXI,iPcls) = pAiv(AIV_PLAG_INDEXI,iPcls)
!    pAivOld(AIV_PLAG_INDEXJ,iPcls) = pAiv(AIV_PLAG_INDEXJ,iPcls)
!    pAivOld(AIV_PLAG_INDEXK,iPcls) = pAiv(AIV_PLAG_INDEXK,iPcls)

    pAiv(AIV_PLAG_ICELLS,iPcls) = indexNew(1)
    pAiv(AIV_PLAG_INDEXI,iPcls) = indexNew(2) 
    pAiv(AIV_PLAG_INDEXJ,iPcls) = indexNew(3) 
    pAiv(AIV_PLAG_INDEXK,iPcls) = indexNew(4)
  
  END DO  ! iPcls 
   
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1010 FORMAT(A,'Time = ', 1PE12.5, ', iReg = ', I5, ', iPcls = ', I8, &
            ', aivOld = ',6I5,', posOld = ',3(1PE15.7),&
            ', posCurr = ', 3(1PE15.7),', pVel = ', 3(1PE15.7), &
            ', mVel = ', 3(1PE15.7),&
            ', diam = ', 1PE15.7,', mass = ',1PE15.7,' , tauL = ', 1PE15.7)

END SUBROUTINE PLAG_getCellIndices
!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_GetCellIndices.F90,v $
! Revision 1.4  2009/10/26 00:19:32  mtcampbe
! Updates for completion of NATIVE_MP_IO
!
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:33  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/08/05 20:10:41  fnajjar
! Included deletion flag for aiv after robust cell search algorithm
!
! Revision 1.8  2004/07/01 14:15:54  fnajjar
! Added variables to IO when particle search fails
!
! Revision 1.7  2004/07/01 14:13:38  fnajjar
! Commented out aivOld loading from aiv since it should be done at initial RK-stage
!
! Revision 1.6  2004/04/09 23:10:20  fnajjar
! Added robust kernel for incell testing
!
! Revision 1.5  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.4  2003/04/18 22:49:04  fnajjar
! Bug fix to move IF statement of cellLocation inside ijk DO-loop
!
! Revision 1.3  2003/01/16 20:43:27  f-najjar
! Include iReg in calling sequence
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:15:43  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








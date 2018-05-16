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
! Purpose: Perform in cell test based on robust search algorithm
!
! Description: none.
!
! Input: 
!   region = data of current region
!   posPlag = particle position vector
!   indexCurr = indices for current cell
!        
! Output: indexNew = new cell index
!          cellLocate = logical variable set to TRUE if test successful
!
! Notes: 
!   1. compute distance between posPlag and face centroids
!       for dummy cells and take minimum distance
!   2. if such distance does not exist, search fails
!
!******************************************************************************
!
! $Id: PLAG_InCellTestRobust.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InCellTestRobust( region, posPlag, indexCurr, &
                                  indexNew,cellLocate )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError  
  USE ModParameters
  
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset, &
                            RFLO_GetDimensPhys, &
                            RFLO_GetDimensPhysNodes

  IMPLICIT NONE

#include "Indexing.h"

! *****************************************************************************
! Declarations
! *****************************************************************************

! =============================================================================  
! Arguments 
! =============================================================================  

  TYPE(t_region) :: region

  INTEGER,      DIMENSION(4), INTENT(IN ) :: indexCurr
  INTEGER,      DIMENSION(4), INTENT(OUT) :: indexNew
  
  LOGICAL,                    INTENT(OUT) :: cellLocate  
  
  REAL(RFREAL), DIMENSION(3), INTENT(IN ) :: posPlag  

! =============================================================================  
! Locals
! =============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: i,iCOff,ijCOff,ijkC,iNOff,ijNOff,ijkNR,ijkNRI,ijkNRJ,ijkNRK,    &
             iLev,ipcbeg,ipcend,ipnbeg,ipnend,j,jpcbeg,jpcend,jpnbeg,jpnend, &
             k,kpcbeg,kpcend,kpnbeg,kpnend,nSearch
  INTEGER,      DIMENSION(4) :: indexSearch
  
  REAL(RFREAL)               :: distPos,distPosMin
  REAL(RFREAL), DIMENSION(3) :: faceCentroid
  
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pFc
  
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InCellTestRobust.F90,v $ $Revision: 1.3 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InCellTestRobust',&
  'PLAG_InCellTestRobust.F90' )

! *****************************************************************************
! Set variables and pointers
! *****************************************************************************
    
  cellLocate = .FALSE.
  distPosMin = HUGE(1.0_RFREAL)
  nSearch    = 0
  indexNew(1:4) = CRAZY_VALUE_INT

  iLev = region%currLevel

  pFc => region%levels(iLev)%plag%fc

! *****************************************************************************
! Get extent of physical cells and nodes
! *****************************************************************************
   
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,  &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,  &
                                jpnbeg,jpnend,kpnbeg,kpnend )

! *****************************************************************************
! Load indices for nodes and cells
! *****************************************************************************
   
  i = indexCurr(2)
  j = indexCurr(3)
  k = indexCurr(4)
   
  ijkNR  = IndIJK(i,  j  ,k  ,iNOff,ijNOff)              
  ijkNRI = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
  ijkNRJ = IndIJK(i,j+1  ,k  ,iNOff,ijNOff)
  ijkNRK = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)   

  ijkC   = IndIJK(i,j,k,iCOff,ijCOff)

! *****************************************************************************
! Perform search for indices of dummy cells
! *****************************************************************************
   
  IF ( i-1 < ipcbeg ) THEN
    faceCentroid(1:3) = pFc(XCOORD:ZCOORD,ICOORD,ijkNR)
    distPos = SUM( (faceCentroid(1:3)-posPlag(1:3))**2 )

    IF ( distPos < distPosMin ) THEN
      distPosMin = distPos
      ijkC = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
      indexSearch(1:4) = (/ijkC,i-1,j,k/)
      nSearch = nSearch+1
    ENDIF ! distPos 

  ENDIF ! i-1
  
  IF ( i+1 > ipcend ) THEN
    faceCentroid(1:3) = pFc(XCOORD:ZCOORD,ICOORD,ijkNRI)
    distPos = SUM( (faceCentroid(1:3)-posPlag(1:3))**2 )

    IF ( distPos < distPosMin ) THEN
      distPosMin = distPos
      ijkC = IndIJK(i+1,j  ,k  ,iCOff,ijCOff)
      indexSearch(1:4) = (/ijkC,i+1,j,k/)
      nSearch = nSearch+1
    ENDIF ! distPos 

  ENDIF ! i+1

  IF ( j-1 < jpcbeg ) THEN
    faceCentroid(1:3) = pFc(XCOORD:ZCOORD,JCOORD,ijkNR)
    distPos = SUM( (faceCentroid(1:3)-posPlag(1:3))**2 )

    IF ( distPos < distPosMin ) THEN
      distPosMin = distPos
      ijkC = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
      indexSearch(1:4) = (/ijkC,i,j-1,k/)
      nSearch = nSearch+1
    ENDIF ! distPos 

  ENDIF ! j-1
  
  IF ( j+1 >  jpcend ) THEN
    faceCentroid(1:3) = pFc(XCOORD:ZCOORD,JCOORD,ijkNRJ)
    distPos = SUM( (faceCentroid(1:3)-posPlag(1:3))**2 )

    IF ( distPos < distPosMin ) THEN
      distPosMin = distPos
      ijkC = IndIJK(i  ,j+1,k  ,iCOff,ijCOff)
      indexSearch(1:4) = (/ijkC,i,j+1,k/)
      nSearch = nSearch+1
    ENDIF ! distPos 

  ENDIF ! j+1   

  IF ( k-1 < kpcbeg ) THEN
    faceCentroid(1:3) = pFc(XCOORD:ZCOORD,KCOORD,ijkNR)
    distPos = SUM( (faceCentroid(1:3)-posPlag(1:3))**2 )

    IF ( distPos < distPosMin ) THEN
      distPosMin = distPos
      ijkC = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
      indexSearch(1:4) = (/ijkC,i,j,k-1/)
      nSearch = nSearch+1
    ENDIF ! distPos 

  ENDIF ! k-1
  
  IF ( k+1 >  kpcend ) THEN
    faceCentroid(1:3) = pFc(XCOORD:ZCOORD,KCOORD,ijkNRK)
    distPos = SUM( (faceCentroid(1:3)-posPlag(1:3))**2 )

    IF ( distPos < distPosMin ) THEN
      distPosMin = distPos
      ijkC = IndIJK(i  ,j  ,k+1,iCOff,ijCOff)
      indexSearch(1:4) = (/ijkC,i,j,k+1/)
      nSearch = nSearch+1
    ENDIF ! distPos 

  ENDIF ! k+1   

  IF ( nSearch /= 0 ) THEN
    cellLocate = .TRUE.
    indexNew(1:4) = indexSearch(1:4)
  ENDIF ! nSearch

! *****************************************************************************
! End
! *****************************************************************************

999  CONTINUE
  CALL DeregisterFunction( global )
 
END SUBROUTINE PLAG_InCellTestRobust

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InCellTestRobust.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/04/09 23:17:19  fnajjar
! Initial Import for robust cell search algorithm
!
!******************************************************************************








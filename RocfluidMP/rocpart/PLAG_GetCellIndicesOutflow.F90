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
! Purpose: Search algorithm for cell indices
!          for outflow boundary conditions
!          since no geometry information is available for dummy cells. 
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%plag%aiv
!         region%levels%plag%aivOld
!
! Notes: Test only if particle is no longer in inner cells.
!
!******************************************************************************
!
! $Id: PLAG_GetCellIndicesOutflow.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_GetCellIndicesOutflow( region )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetNodeOffset,   RFLO_GetCellOffset,     &
                            RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, iLev, lbound, nPatches, nPcls
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkNPatch
  INTEGER :: iCPlag, ipcbeg, ipcend, ibeg, iend, idir
  INTEGER :: jCPlag, jpcbeg, jpcend, jbeg, jend, jdir
  INTEGER :: kCPlag, kpcbeg, kpcend, kbeg, kend, kdir

  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

  LOGICAL :: lboundSkip(6)

  REAL(RFREAL)               :: dpGrid, sgn
  REAL(RFREAL), DIMENSION(3) :: diffPos, faceCentroid, posPlag, sFace
  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: pCv, pSi, pSj, pSk
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pFc

  TYPE(t_patch),  POINTER :: pPatch
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_GetCellIndicesOutflow.F90,v $ $Revision: 1.3 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_GetCellIndicesOutflow',&
  'PLAG_GetCellIndicesOutflow.F90' )

! Get dimensions --------------------------------------------------------------

  iLev  = region%currLevel
  nPcls = region%levels(iLev)%plag%nPcls 
  
  IF (nPcls == 0) GOTO 999 ! exit if no particles

  nPatches = region%nPatches

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! Set pointers ----------------------------------------------------------------

  pPlag   => region%levels(iLev)%plag 
  pAiv    => pPlag%aiv 
  pAivOld => pPlag%aivOld 
  pCv     => pPlag%cv

  pSi     => region%levels(iLev)%plag%si
  pSj     => region%levels(iLev)%plag%sj
  pSk     => region%levels(iLev)%plag%sk

  pFc     => region%levels(iLev)%plag%fc 

! Loop over Lagrangian particles ----------------------------------------------

  DO iPcls = 1, nPcls

    iCPlag = pAiv(AIV_PLAG_INDEXI,iPcls)
    jCPlag = pAiv(AIV_PLAG_INDEXJ,iPcls)
    kCPlag = pAiv(AIV_PLAG_INDEXK,iPcls)
    
! - Cycle if particle is not adjacent to a boundary ---------------------------

    IF (ipcbeg < iCPlag .AND. iCPlag < ipcend .AND. &
        jpcbeg < jCPlag .AND. jCPlag < jpcend .AND. &
        kpcbeg < kCPlag .AND. kCPlag < kpcend ) CYCLE

    posPlag(XCOORD:ZCOORD) = pCv(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls)
        
! - Set lboundSkip(:) to .TRUE. for values of lbound adjacent to a boundary ---

    lboundSkip(1:6) = .TRUE.

    IF (iCPlag <= ipcbeg) lboundSkip(1) = .FALSE. ! to include ipcbeg = ipcend
    IF (iCPlag >= ipcend) lboundSkip(2) = .FALSE. ! case, ELSE IF not used here

    IF (jCPlag <= jpcbeg) lboundSkip(3) = .FALSE.
    IF (jCPlag >= jpcend) lboundSkip(4) = .FALSE.

    IF (kCPlag <= kpcbeg) lboundSkip(5) = .FALSE.
    IF (kCPlag >= kpcend) lboundSkip(6) = .FALSE.

! - Loop over patches ---------------------------------------------------------

    DO iPatch=1,nPatches

      pPatch => region%levels(iLev)%patches(iPatch)
      lbound =  pPatch%lbound
      IF (lboundSkip(lbound)) CYCLE

      bcType = pPatch%bcType
    
! - Select outflow boundary condition type ------------------------------------

      IF ( bcType>=BC_OUTFLOW .AND. bcType<=BC_OUTFLOW+BC_RANGE ) THEN 

! ----- Check if particle cell is within (physical cells of) patch ------------

        CALL RFLO_GetPatchIndices( region,pPatch,iLev, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )

        IF ( ibeg <= iCPlag .AND. iCPlag <= iend .AND. &
             jbeg <= jCPlag .AND. jCPlag <= jend .AND. &
             kbeg <= kCPlag .AND. kCPlag <= kend ) THEN

! ------- Select correct face vector and make it point inwards ----------------

          CALL RFLO_GetPatchDirection( pPatch,idir,jdir,kdir )

          IF (lbound==1 .OR. lbound==3 .OR. lbound==5) THEN
            sgn   = -1.0_RFREAL
            ijkNPatch = IndIJK(iCPlag,jCPlag,kCPlag,iNOff,ijNOff)
          ELSE
            sgn   = +1.0_RFREAL
            ijkNPatch = IndIJK(iCPlag-idir,jCPlag-jdir,kCPlag-kdir,iNOff,ijNOff)
          ENDIF ! lbound

          SELECT CASE (lbound)

            CASE(1,2)
              sFace(1:3)        = sgn*pSi(XCOORD:ZCOORD,       ijkNPatch)
              faceCentroid(1:3) =     pFc(XCOORD:ZCOORD,ICOORD,ijkNPatch)

            CASE(3,4)
              sFace(1:3)        = sgn*pSj(XCOORD:ZCOORD,       ijkNPatch)
              faceCentroid(1:3) =     pFc(XCOORD:ZCOORD,JCOORD,ijkNPatch)

            CASE(5,6)
              sFace(1:3)        = sgn*pSk(XCOORD:ZCOORD,       ijkNPatch)
              faceCentroid(1:3) =     pFc(XCOORD:ZCOORD,KCOORD,ijkNPatch)

          END SELECT ! lbound

          diffPos(1:3) = posPlag(1:3)-faceCentroid(1:3)
          dpGrid       = DOT_PRODUCT(sFace,diffPos)

! ------- Test for particle exiting computational domain ----------------------

          IF ( dpGrid < 0.0_RFREAL ) THEN

            pAivOld(AIV_PLAG_ICELLS,iPcls) = pAiv(AIV_PLAG_ICELLS,iPcls)
            pAivOld(AIV_PLAG_INDEXI,iPcls) = iCPlag
            pAivOld(AIV_PLAG_INDEXJ,iPcls) = jCPlag
            pAivOld(AIV_PLAG_INDEXK,iPcls) = kCPlag

            iCPlag = iCPlag - idir ! update indices to new (dummy) cell
            jCPlag = jCPlag - jdir
            kCPlag = kCPlag - kdir

            pAiv(   AIV_PLAG_ICELLS,iPcls) = IndIJK(iCPlag,jCPlag,kCPlag,iCOff,ijCOff)
            pAiv(   AIV_PLAG_INDEXI,iPcls) = iCPlag
            pAiv(   AIV_PLAG_INDEXJ,iPcls) = jCPlag
            pAiv(   AIV_PLAG_INDEXK,iPcls) = kCPlag

          ENDIF ! dpGrid
        ENDIF   ! iCPlag, jCPlag, kCPlag
      ENDIF     ! bcType
    ENDDO       ! iPatch
  ENDDO         ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_GetCellIndicesOutflow

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_GetCellIndicesOutflow.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:34  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/11/03 21:21:51  fnajjar
! Changed definition of face vectors pointing to PLAG datastructure
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/05/01 22:58:30  jferry
! overhauled structure in order to optimize performance
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:15:43  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************








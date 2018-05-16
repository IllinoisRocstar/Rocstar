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
! Purpose: remove particle from datastructure if exiting computational domain. 
!
! Description: none.
!
! Input: region = current region
!        iReg   = number of current region.
!
! Output: regions(iReg)%levels%patch%buffPlag%aiv,arv,cv,dv,tv = buffer data.
!         regions(iReg)%levels%plag%aiv,arv,cv,dv,tv           = Plag data.
!
! Notes:
!
!   Particle indices have been updated at the outflow BC only, so particles
!   are removed if their indices are no longer within the region
!
!******************************************************************************
!
! $Id: PLAG_PatchRemoveDataOutflow.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PatchRemoveDataOutflow( region, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER        :: iReg

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nPcls, nPclsPrev, iPclsRegionIn
  INTEGER :: iCPlag, ipcbeg, ipcend
  INTEGER :: jCPlag, jpcbeg, jpcend
  INTEGER :: kCPlag, kpcbeg, kpcend

  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

  LOGICAL :: lboundSkip(6)

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvOld, pCv, pCvOld, &
                                           pDv, pTv, pRhs, pRhsSum

  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_PatchRemoveDataOutflow.F90,v $ $Revision: 1.3 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_PatchRemoveDataOutflow',&
  'PLAG_PatchRemoveDataOutflow.F90' )

! Get dimensions --------------------------------------------------------------

  iLev  = region%currLevel
  nPcls = region%levels(iLev)%plag%nPcls 

  IF (nPcls == 0) GOTO 999 ! exit if no particles

! Set pointers ----------------------------------------------------------------

  pPlag   => region%levels(iLev)%plag   
  pAiv    => pPlag%aiv
  pArv    => pPlag%arv
  pCv     => pPlag%cv
  pDv     => pPlag%dv
  pTv     => pPlag%tv
  pRhs    => pPlag%rhs
  pRhsSum => pPlag%rhsSum

  pAivOld => pPlag%aivOld
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld

! Get grid dimensions ---------------------------------------------------------
  
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
     
! Initialize counters for particles inside and outside region -----------------

  iPclsRegionIn = 0

! Loop over particles ---------------------------------------------------------

  DO iPcls = 1, nPcls

    iCPlag = pAiv(AIV_PLAG_INDEXI,iPcls)
    jCPlag = pAiv(AIV_PLAG_INDEXJ,iPcls)
    kCPlag = pAiv(AIV_PLAG_INDEXK,iPcls)
    
! - Check if particle cell is within region -----------------------------------

    IF ( ipcbeg <= iCPlag .AND. iCPlag <= ipcend .AND. &
         jpcbeg <= jCPlag .AND. jCPlag <= jpcend .AND. &
         kpcbeg <= kCPlag .AND. kCPlag <= kpcend ) THEN

! --- Shift particle datastructure only if particle is not in its proper spot -
 
      iPclsRegionIn = iPclsRegionIn + 1 

      IF ( iPclsRegionIn /= iPcls ) THEN
        pAiv(   :,iPclsRegionIn) = pAiv(   :,iPcls)
        pArv(   :,iPclsRegionIn) = pArv(   :,iPcls)
        pCv(    :,iPclsRegionIn) = pCv(    :,iPcls)
        pDv(    :,iPclsRegionIn) = pDv(    :,iPcls)
        pTv(    :,iPclsRegionIn) = pTv(    :,iPcls)
        pRhs(   :,iPclsRegionIn) = pRhs(   :,iPcls)
        pRhsSum(:,iPclsRegionIn) = pRhsSum(:,iPcls)

        pAivOld(:,iPclsRegionIn) = pAivOld(:,iPcls)
        pArvOld(:,iPclsRegionIn) = pArvOld(:,iPcls)
        pCvOld( :,iPclsRegionIn) = pCvOld( :,iPcls)
      ENDIF ! iPclsRegionIn
 
    ENDIF   ! iCPlag, jCPlag, kCPlag

  ENDDO     ! iPcls

! Get new particle datasize ---------------------------------------------------

  nPclsPrev   = pPlag%nPcls    
  pPlag%nPcls = iPclsRegionIn

!  WRITE(STDOUT,'(A,I4,2I8)') &
!  '      PLAG_PatchRemoveDataOutflow: iReg, nPclsPrev,nPclsCurr',iReg,nPclsPrev,pPlag%nPcls
      
! finalize --------------------------------------------------------------------

999  CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_PatchRemoveDataOutflow

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchRemoveDataOutflow.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:58  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2003/05/27 19:19:57  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.7  2003/05/01 22:58:30  jferry
! overhauled structure in order to optimize performance
!
! Revision 1.6  2003/05/01 18:07:25  jferry
! intermediate form checked in for documentation purposes only
!
! Revision 1.5  2003/04/17 00:00:27  fnajjar
! Commented out I/O to STDOUT
!
! Revision 1.4  2003/01/17 20:12:32  f-najjar
! Fixed FORMAT statements with iReg
!
! Revision 1.3  2003/01/17 19:30:24  f-najjar
! Added iReg in calling sequence
!
! Revision 1.2  2003/01/16 20:24:32  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:18:35  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************








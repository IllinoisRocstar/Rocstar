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
! Purpose: append Lagrangian particle datastructure from buffers 
!          communicated via the patches. 
!
! Description: none.
!
! Input: regions = regions
!        iReg    = current region number.
!
! Output: regions(iReg)%levels%plag%aiv,arv,cv,dv,tv = Plag data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_AppendDataFromBuffers.F90,v 1.4 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_AppendDataFromBuffers( region, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection
  
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER        :: iReg

! ... loop variables
  INTEGER :: iPatch, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, ibeg, iBuff, idir, iend, iLev,               &
             jbeg, jdir, jend, kbeg, kdir, kend,                  &
             lbound, n1, n2, nBuffSizeDes, nOff, nPatches,        &
             nPcls, nPclsBuffTot, nPclsEnd, nPclsPrev, nPclsStart 
   
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv,  pAivBuff, pAivOld, pAivOldBuff

  LOGICAL :: plagRegionIn, plagPatchDumCell

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvBuff,              &
                                           pArvOld, pArvOldBuff,        &
                                           pCv, pCvBuff,                &
                                           pCvOld, pCvOldBuff,          &
                                           pDv, pDvBuff,                &
                                           pRhs, pRhsBuff, pRhsSum,     &
                                           pRhsSumBuff, pTv, pTvBuff
  
  TYPE(t_patch),       POINTER :: pPatch
  TYPE(t_plag),        POINTER :: pPlag 
  TYPE(t_buffer_plag), POINTER :: pBuffPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_AppendDataFromBuffers.F90,v $ $Revision: 1.4 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_AppendDataFromBuffers',&
  'PLAG_AppendDataFromBuffers.F90' )

! Get dimensions --------------------------------------------------------------

  iLev     = region%currLevel
  nPatches = region%nPatches

  nPcls        = region%levels(iLev)%plag%nPcls 
  nPclsBuffTot = region%plagInput%nPclsBuffTot
  nPclsPrev    = 0
  
! Set pointers ----------------------------------------------------------------
   
  pPlag   => region%levels(iLev)%plag   
  pAiv    => pPlag%aiv
  pArv    => pPlag%arv
  pCv     => pPlag%cv
  pDv     => pPlag%dv
  pTv     => pPlag%tv
  pRhs    => pPlag%rhs
  pRhsSum => pPlag%rhsSum

  pAivOld    => pPlag%aivOld
  pArvOld    => pPlag%arvOld
  pCvOld     => pPlag%cvOld

! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    pPatch  => region%levels(iLev)%patches(iPatch)
    bcType  = pPatch%bcType    
   
! - Get patch dimensions -----------------------------------------------------

    lbound = pPatch%lbound

    CALL RFLO_GetPatchIndices( region,pPatch,iLev, &
                               ibeg,iend,jbeg,jend,kbeg,kend )
    CALL RFLO_GetPatchDirection( pPatch,idir,jdir,kdir )
       
! - Select communication boundary conditions ---------------------------------
        
    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
        (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN 

        pBuffPlag => pPatch%bufferPlag 
        pAivBuff  => pBuffPlag%aiv
        pArvBuff  => pBuffPlag%arv
        pCvBuff   => pBuffPlag%cv
        pDvBuff   => pBuffPlag%dv
        pTvBuff   => pBuffPlag%tv
        pRhsBuff  => pBuffPlag%rhs
        pRhsSumBuff  => pBuffPlag%rhsSum

        pAivOldBuff  => pBuffPlag%aivOld
        pArvOldBuff  => pBuffPlag%arvOld
        pCvOldBuff   => pBuffPlag%cvOld
    
        nPclsStart = 0
        nPclsEnd   = 0
        nPclsPrev  = nPcls
        
! - Get buffer size and start appending for non-null size ---------------------

        nBuffSizeDes = pBuffPlag%nBuffSizeDes

        IF ( nBuffSizeDes /= 0 ) THEN 
        
! - Set loop extent-------------------------------------------------------------

          nPclsStart = nPcls+1
          nPclsEnd   = nPclsStart + (nBuffSizeDes-1)
          
#ifdef PLAG_DEBUG
          WRITE(STDOUT,'(A,I2,3I6)') &
          '      PLAG_AppendDataFromBuffers-iReg: nBuffSizeDes,nPclsStart,nPclsEnd',&
           iReg, nBuffSizeDes,nPclsStart,nPclsEnd
#endif


          DO iPcls = nPclsStart, nPclsEnd
          
            iBuff = iPcls-nPclsStart+1
 
            pAiv(:,iPcls) = pAivBuff(:,iBuff)
            pArv(:,iPcls) = pArvBuff(:,iBuff)
            pCv(:,iPcls)  = pCvBuff( :,iBuff)
            pDv(:,iPcls)  = pDvBuff( :,iBuff)
            pTv(:,iPcls)  = pTvBuff( :,iBuff)
            pRhs(:,iPcls) = pRhsBuff(:,iBuff)
            pRhsSum(:,iPcls) = pRhsSumBuff(:,iBuff)

            pAivOld(:,iPcls) = pAivOldBuff(:,iBuff)
            pArvOld(:,iPcls) = pArvOldBuff(:,iBuff)
            pCvOld(:,iPcls)  = pCvOldBuff( :,iBuff)
            
          END DO ! iPcls

! - Get new particle datasize and reset loop extent----------------------------
          
          nPcls = nPcls+nBuffSizeDes
          region%levels(iLev)%plag%nPcls = nPcls

        END IF ! nBuffSizeDes

    END IF ! bcType 

  ENDDO ! iPatch

#ifdef PLAG_DEBUG
  IF( nPclsPrev /= pPlag%nPcls ) &
  WRITE(STDOUT,'(A,I4,2I8)') &
  '      PLAG_AppendDataFromBuffers-iReg: nPclsPrev,nPclsCurr',&
  iReg,nPclsPrev,pPlag%nPcls
#endif 

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_AppendDataFromBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_AppendDataFromBuffers.F90,v $
! Revision 1.4  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:56:53  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/04/09 23:01:42  fnajjar
! Put WRITE statements inside ifdef construct
!
! Revision 1.8  2003/05/27 19:14:16  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.7  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.6  2003/05/07 00:02:34  fnajjar
! Commented out I/O for nPcls
!
! Revision 1.5  2003/05/06 23:42:46  fnajjar
! Removed exit trap for zero nPcls to allow particles to be communicated
!
! Revision 1.4  2003/04/17 00:00:27  fnajjar
! Commented out I/O to STDOUT
!
! Revision 1.3  2003/02/28 21:44:13  fnajjar
! Initialize nPclsPrev and Exit for Null nPcls
!
! Revision 1.2  2003/01/23 17:21:57  f-najjar
! Removed Hidden TABs
!
! Revision 1.1  2003/01/16 22:42:09  f-najjar
! Initial import
!
!******************************************************************************








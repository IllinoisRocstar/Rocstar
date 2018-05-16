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
! Purpose: receive buffer data from adjacent regions on different processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: integer and real data buffer from other processors.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_BufferDataRecv.F90,v 1.6 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_BufferDataRecv( regions, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  
  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch, iBuff, iCont

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef MPI
  INTEGER :: statusPlag(MPI_STATUS_SIZE)
#endif

  INTEGER :: bcType, iLev, iPatchDes, iPatchSrc, iRegSrc,   &
             nAiv, nArv, nBuffSizeDes, nCont, nCv, nDimI,   &
             nDimR, nDv, nPatches, nRecvBuffI, nRecvBuffR,  &
             nTv, procDes, procSrc, tagSrcI, tagSrcR
  
  INTEGER :: iAiv, iArv, iArvOld, iBuffRecv, iCv, iCvOld,   &
             iRhs, iRhsSum
    
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pRecvBuffI 

  INTEGER, POINTER, DIMENSION(:,:) :: pAivDes, pAivOldDes

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pRecvBuffR
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvDes, pArvOldDes, &
                                           pCvDes , pCvOldDes , &
                                           pRhsDes, pRhsSumDes
                                   
  TYPE(t_patch),  POINTER :: pPatchSrc, pPatchDes
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_BufferDataRecv.F90,v $ $Revision: 1.6 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_BufferDataRecv',&
  'PLAG_BufferDataRecv.F90' )
            
! receive data buffer from source to destination region 
  
! get dimensions --------------------------------------------------------------

  nCont = regions(iReg)%plagInput%nCont
  nCv   = CV_PLAG_LAST + nCont
  nDv   = DV_PLAG_LAST
  nTv   = TV_PLAG_LAST
  nAiv  = AIV_PLAG_LAST
  nArv  = ARV_PLAG_LAST
       
  nDimI = nAiv
  nDimR = 2*nArv +4*nCv     

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches
  
  pPlag       => regions(iReg)%levels(iLev)%plag
  pCvPlagMass => pPlag%cvPlagMass

! loop over patches -----------------------------------------------------------

  DO iPatch = 1, nPatches

! - pointer is at Des region getting data from Src region ---------------------
     
    pPatchDes => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = pPatchDes%bcType
    iRegSrc   = pPatchDes%srcRegion
    iPatchSrc = pPatchDes%srcPatch

! - region interface for various boundary conditions --------------------------

    IF ( (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
         (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
         (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
         (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
         (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) ) THEN

      IF ( regions(iRegSrc)%procid /= global%myProcid ) THEN
        pPatchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
        
        procDes = global%myProcid

! -- set pointers and buffer size on destination region ---------------------------
   
        pAivDes    => pPatchDes%bufferPlag%aiv
        pArvDes    => pPatchDes%bufferPlag%arv

        pCvDes     => pPatchDes%bufferPlag%cv
        pRhsDes    => pPatchDes%bufferPlag%rhs
        pRhsSumDes => pPatchDes%bufferPlag%rhsSum
   
        pAivOldDes => pPatchDes%bufferPlag%aivOld
        pArvOldDes => pPatchDes%bufferPlag%arvOld
        pCvOldDes  => pPatchDes%bufferPlag%cvOld 
                   
        pRecvBuffI => pPatchDes%bufferPlag%recvBuffI
        pRecvBuffR => pPatchDes%bufferPlag%recvBuffR        

        nBuffSizeDes = pPatchDes%bufferPlag%nBuffSizeDes 

! -- exit for null buffer size ------------------------------------------------
        
        IF ( nBuffSizeDes == 0 ) GOTO 999

! -- set buffer sizes ---------------------------------------------------------

        nRecvBuffI = nDimI * nBuffSizeDes
        nRecvBuffR = nDimR * nBuffSizeDes

! -- receive data -------------------------------------------------------------
    
#ifdef MPI

        procSrc = regions(iRegSrc)%procid
        
! --- integer -----------------------------------------------------------------

        tagSrcI = regions(iReg)%localNumber &
                + PLAG_TAG_SHIFT +MPI_PATCHOFF*pPatchSrc%srcPatch*iReg + procDes +1
        
        IF(tagSrcI .gt. global%mpiTagMax) tagSrcI = MOD(tagSrcI,global%mpiTagMax)
        CALL MPI_Recv( pRecvBuffI, nRecvBuffI, MPI_INTEGER,  &
                       procSrc, tagSrcI, global%mpiComm,     &
                       statusPlag, global%mpierr )

        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

#ifdef PLAG_MPI_DEBUG                    
   IF(nRecvBuffI /=0 ) &
   WRITE(STDOUT,*) '  PLAG_BufferDataRecv-INT: iRegDes, iRegSrc, procSrc, tagSrcI, nRecvBuffI  = ',&
                      iReg, iRegSrc, procSrc,tagSrcI, nRecvBuffI
#endif

! --- real --------------------------------------------------------------------

        tagSrcR = regions(iReg)%localNumber &
                + PLAG_TAG_SHIFT +MPI_PATCHOFF*pPatchSrc%srcPatch*iReg + procDes +2
        
        IF(tagSrcR .gt. global%mpiTagMax) tagSrcR = MOD(tagSrcR,global%mpiTagMax)
        CALL MPI_Recv( pRecvBuffR, nRecvBuffR, MPI_RFREAL, &
                       procSrc, tagSrcR, global%mpiComm,   &
                       statusPlag,global%mpierr )

        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

#ifdef PLAG_MPI_DEBUG                    
   IF(nRecvBuffR /=0 ) &
   WRITE(STDOUT,*) '  PLAG_BufferDataRecv-REAL: iRegDes, iRegSrc, procSrc, tagSrcR, nRecvBuffR  = ',&
                  iReg, iRegSrc, procSrc,tagSrcR, nRecvBuffR
#endif

#endif

! -- copy receive buffers into local arrays -----------------------------------

! --- integer variables -------------------------------------------------------
      
        DO iBuff = 1, nBuffSizeDes 

! ---- compute shift ----------------------------------------------------------

          iBuffRecv = nDimI*(iBuff-1) +1
          iAiv = iBuffRecv

          pAivDes(AIV_PLAG_PIDINI,iBuff) = pRecvBuffI(iAiv  )
          pAivDes(AIV_PLAG_REGINI,iBuff) = pRecvBuffI(iAiv+1) 
          pAivDes(AIV_PLAG_REGCRT,iBuff) = pRecvBuffI(iAiv+2) 
          pAivDes(AIV_PLAG_ICELLS,iBuff) = pRecvBuffI(iAiv+3) 
          pAivDes(AIV_PLAG_INDEXI,iBuff) = pRecvBuffI(iAiv+4)
          pAivDes(AIV_PLAG_INDEXJ,iBuff) = pRecvBuffI(iAiv+5)
          pAivDes(AIV_PLAG_INDEXK,iBuff) = pRecvBuffI(iAiv+6)
          pAivDes(AIV_PLAG_BURNSTAT,iBuff) = pRecvBuffI(iAiv+7) 
          pAivDes(AIV_PLAG_STATUS,iBuff) = pRecvBuffI(iAiv+8) 
          
          pAivOldDes(AIV_PLAG_PIDINI,iBuff) = pRecvBuffI(iAiv  )
          pAivOldDes(AIV_PLAG_REGINI,iBuff) = pRecvBuffI(iAiv+1) 
          pAivOldDes(AIV_PLAG_REGCRT,iBuff) = pRecvBuffI(iAiv+2) 
          pAivOldDes(AIV_PLAG_ICELLS,iBuff) = pRecvBuffI(iAiv+3) 
          pAivOldDes(AIV_PLAG_INDEXI,iBuff) = pRecvBuffI(iAiv+4)
          pAivOldDes(AIV_PLAG_INDEXJ,iBuff) = pRecvBuffI(iAiv+5)
          pAivOldDes(AIV_PLAG_INDEXK,iBuff) = pRecvBuffI(iAiv+6)
          pAivOldDes(AIV_PLAG_BURNSTAT,iBuff) = pRecvBuffI(iAiv+7)
          pAivOldDes(AIV_PLAG_STATUS,iBuff) = pRecvBuffI(iAiv+8)

#ifdef PLAG_MPI_DEBUG                   
   IF(nBuffSizeDes /=0 ) &
  WRITE(STDOUT,*) '  PLAG_BufferDataRecv-INT: procDes, iBuff, iAiv, pAivDes = ',&
                     procDes, iBuff, iAiv, pAivDes(:,iBuff)
#endif
 
        ENDDO ! iBuff

! --- real variables ----------------------------------------------------------
      
        DO iBuff = 1, nBuffSizeDes 

! ---- compute shifts ---------------------------------------------------------

          iBuffRecv = nDimR*(iBuff-1) +1
          iCv       = iBuffRecv
          iRhs      = iBuffRecv +nCv 
          iRhsSum   = iBuffRecv +2*nCv 
          iCvOld    = iBuffRecv +3*nCv 
          iArv      = iBuffRecv +4*nCv 
          iArvOld   = iBuffRecv +4*nCv +nArv

#ifdef PLAG_MPI_DEBUG 
   IF(nBuffSizeDes /=0 ) &
  WRITE(STDOUT,*) '  PLAG_BufferDataRecv-REAL: procDes,iBuff, iBuffRecv, iCv, iRhs, iRhsSum, iCvOld, iArv, iArvOld = ',&
                     procDes, iBuff, iBuffRecv,iCv, iRhs, iRhsSum, iCvOld, iArv, iArvOld  
  
   IF(nBuffSizeDes /=0 ) &
  WRITE(STDOUT,*) '  PLAG_BufferDataRecv-Entering CV:procDes, iReg, iBuff, pRecvBuffR(iCv)', &
                     procDes,iReg, iBuff,pRecvBuffR(iCv:iCv+nCv)  
#endif
            
! ---- load cv ----------------------------------------------------------------

          pCvDes(CV_PLAG_XMOM,iBuff) = pRecvBuffR(iCv  )  
          pCvDes(CV_PLAG_YMOM,iBuff) = pRecvBuffR(iCv+1)  
          pCvDes(CV_PLAG_ZMOM,iBuff) = pRecvBuffR(iCv+2)  
          pCvDes(CV_PLAG_ENER,iBuff) = pRecvBuffR(iCv+3) 
          pCvDes(CV_PLAG_XPOS,iBuff) = pRecvBuffR(iCv+4)  
          pCvDes(CV_PLAG_YPOS,iBuff) = pRecvBuffR(iCv+5)  
          pCvDes(CV_PLAG_ZPOS,iBuff) = pRecvBuffR(iCv+6)
          pCvDes(CV_PLAG_ENERVAPOR,iBuff) = pRecvBuffR(iCv+7)   
          DO iCont = 1, nCont
            pCvDes(pCvPlagMass(iCont),iBuff) = pRecvBuffR(iCv+(CV_PLAG_LAST-1)+iCont)  
          ENDDO ! iCont

#ifdef PLAG_MPI_DEBUG        
   IF(nBuffSizeDes /=0 ) &
         WRITE(STDOUT,*) '  PLAG_BufferDataRecv-Done with CV:procDes, iReg, iBuff, pCvDes', &
                          procDes,iReg, iBuff,pCvDes(:,iBuff)
#endif
          
! ---- load rhs ---------------------------------------------------------------          


          pRhsDes(CV_PLAG_XMOM,iBuff) = pRecvBuffR(iRhs  ) 
          pRhsDes(CV_PLAG_YMOM,iBuff) = pRecvBuffR(iRhs+1)   
          pRhsDes(CV_PLAG_ZMOM,iBuff) = pRecvBuffR(iRhs+2)   
          pRhsDes(CV_PLAG_ENER,iBuff) = pRecvBuffR(iRhs+3)   
          pRhsDes(CV_PLAG_XPOS,iBuff) = pRecvBuffR(iRhs+4)  
          pRhsDes(CV_PLAG_YPOS,iBuff) = pRecvBuffR(iRhs+5)  
          pRhsDes(CV_PLAG_ZPOS,iBuff) = pRecvBuffR(iRhs+6) 
          pRhsDes(CV_PLAG_ENERVAPOR,iBuff) = pRecvBuffR(iRhs+7) 
          DO iCont = 1, nCont
            pRhsDes(pCvPlagMass(iCont),iBuff) = pRecvBuffR(iRhs+(CV_PLAG_LAST-1)+iCont)  
          ENDDO ! iCont

#ifdef PLAG_MPI_DEBUG          
   IF(nBuffSizeDes /=0 ) &
         WRITE(STDOUT,*) '  PLAG_BufferDataRecv-Done with RhsSum:procDes, iReg, iBuff, pRhsDes', &
                          procDes,iReg, iBuff, pRhsDes(:,iBuff)
#endif
          
! ---- load rhsSum ------------------------------------------------------------

          
          pRhsSumDes(CV_PLAG_XMOM,iBuff) = pRecvBuffR(iRhsSum  )   
          pRhsSumDes(CV_PLAG_YMOM,iBuff) = pRecvBuffR(iRhsSum+1)   
          pRhsSumDes(CV_PLAG_ZMOM,iBuff) = pRecvBuffR(iRhsSum+2)  
          pRhsSumDes(CV_PLAG_ENER,iBuff) = pRecvBuffR(iRhsSum+3)  
          pRhsSumDes(CV_PLAG_XPOS,iBuff) = pRecvBuffR(iRhsSum+4)  
          pRhsSumDes(CV_PLAG_YPOS,iBuff) = pRecvBuffR(iRhsSum+5)  
          pRhsSumDes(CV_PLAG_ZPOS,iBuff) = pRecvBuffR(iRhsSum+6)
          pRhsSumDes(CV_PLAG_ENERVAPOR,iBuff) = pRecvBuffR(iRhsSum+7)
          DO iCont = 1, nCont
            pRhsSumDes(pCvPlagMass(iCont),iBuff) = pRecvBuffR(iRhsSum+(CV_PLAG_LAST-1)+iCont) 
          ENDDO ! iCont          

#ifdef PLAG_MPI_DEBUG          
   IF(nBuffSizeDes /=0 ) &
         WRITE(STDOUT,*) '  PLAG_BufferDataRecv-Done with RhsSum:procDes, iReg, iBuff, pRhsSumDes', &
                          procDes,iReg, iBuff, pRhsSumDes(:,iBuff)
#endif
          
! ---- load cvOld -------------------------------------------------------------

          
          pCvOldDes(CV_PLAG_XMOM,iBuff) = pRecvBuffR(iCvOld  )  
          pCvOldDes(CV_PLAG_YMOM,iBuff) = pRecvBuffR(iCvOld+1)   
          pCvOldDes(CV_PLAG_ZMOM,iBuff) = pRecvBuffR(iCvOld+2) 
          pCvOldDes(CV_PLAG_ENER,iBuff) = pRecvBuffR(iCvOld+3) 
          pCvOldDes(CV_PLAG_XPOS,iBuff) = pRecvBuffR(iCvOld+4)  
          pCvOldDes(CV_PLAG_YPOS,iBuff) = pRecvBuffR(iCvOld+5) 
          pCvOldDes(CV_PLAG_ZPOS,iBuff) = pRecvBuffR(iCvOld+6)
          pCvOldDes(CV_PLAG_ENERVAPOR,iBuff) = pRecvBuffR(iCvOld+7) 
          DO iCont = 1, nCont
            pCvOldDes(pCvPlagMass(iCont),iBuff) = pRecvBuffR(iCvOld+(CV_PLAG_LAST-1)+iCont) 
          ENDDO ! iCont

#ifdef PLAG_MPI_DEBUG          
   IF(nBuffSizeDes /=0 ) &
         WRITE(STDOUT,*) '  PLAG_BufferDataRecv-Done with CVOld:procDes, iReg, iBuff, pCvDes', &
                          procDes,iReg, iBuff,pCvOldDes(:,iBuff)
#endif
          
! ---- load arv and arvOld ----------------------------------------------------
         
             pArvDes(ARV_PLAG_SPLOAD,iBuff) = pRecvBuffR(iArv  )
             pArvDes(ARV_PLAG_DISTOT,iBuff) = pRecvBuffR(iArv+1)
 
          pArvOldDes(ARV_PLAG_SPLOAD,iBuff) = pRecvBuffR(iArvOld  ) 
          pArvOldDes(ARV_PLAG_DISTOT,iBuff) = pRecvBuffR(iArvOld+1) 

#ifdef PLAG_MPI_DEBUG          
   IF(nBuffSizeDes /=0 ) &
         WRITE(STDOUT,*) '  PLAG_BufferDataRecv-Done with Arv:procDes, iReg, iBuff, pArvDes', &
                          procDes,iReg, iBuff,pArvDes(:,iBuff),pArvOldDes(:,iBuff)
#endif

        ENDDO ! iBuff 
 
      ENDIF ! regions
    ENDIF ! bcType

999 CONTINUE
         
  ENDDO ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_BufferDataRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_BufferDataRecv.F90,v $
! Revision 1.6  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.5  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.2  2005/05/31 21:37:31  fnajjar
! Added ARV_PLAG_DISTOT for proper IO capabilities
!
! Revision 1.1  2004/12/01 20:56:55  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/04/09 23:04:12  fnajjar
! Added AIV_PLAG_STATUS to buffers being sent and received
!
! Revision 1.6  2004/03/21 00:43:32  fnajjar
! Fixed tags to be smaller number since Frost run-time system complains about size
!
! Revision 1.5  2004/03/12 23:41:46  fnajjar
! Bug fix for pRhsSumDes with incorrect assignment
!
! Revision 1.4  2004/03/06 21:25:05  fnajjar
! Added PLAG_TAG_SHIFT to MPI-based communication tags
!
! Revision 1.3  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.2  2003/05/13 15:02:13  fnajjar
! Added IFDEF clause around WRITE statement
!
! Revision 1.1  2003/02/21 17:08:48  fnajjar
! Initial import
!
!******************************************************************************








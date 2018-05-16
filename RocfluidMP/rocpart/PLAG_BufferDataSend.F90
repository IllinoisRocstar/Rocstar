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
! Purpose: send buffer data to adjacent region on different processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: integer and real data buffer to other processors.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_BufferDataSend.F90,v 1.6 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_BufferDataSend( regions, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModIndexing, ONLY   : GetIJK, IndIJKMap
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetPatchMapping
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

  INTEGER :: bcType, iLev, iPatchDes, iRegDes,     &
             iRequestPlag,  nArv, nAiv, nBuffSizeSrc, nCont,  &
             nCv, nDimI, nDimR, nDv, nPatches, nSendBuffI,    &
             nSendBuffR,  nTv, procDes, tagDesI, tagDesR
 
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc,    &
             idirSrc, jdirSrc, kdirSrc, iCOffSrc, ijCOffSrc,          &
             ibegDes, iendDes, jbegDes, jendDes, kbegDes, kendDes,    &
             idirDes, jdirDes, kdirDes, iCOffDes, ijCOffDes, iLevDes
  
  INTEGER :: lbSrc, lbDes, l1DesDir, l2DesDir, mapMat(3,4),           &
             nDumCellsSrc, nDumCellsDes

  INTEGER :: indexISrc, indexJSrc, indexKSrc, iCellsSrc, &
             indexIDes, indexJDes, indexKDes, iCellsDes

  INTEGER :: iAiv, iArv, iArvOld, iBuffSend, iCv, iCvOld, iRhs, iRhsSum
    
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pSendBuffI 

  INTEGER, POINTER, DIMENSION(:,:) :: pAivSrc

  LOGICAL :: alignSrc, alignDes
    
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSendBuffR
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvSrc, pArvOldSrc,  &
                                           pCvSrc , pCvOldSrc ,  &
                                           pRhsSrc, pRhsSumSrc
                                       
  TYPE(t_patch),  POINTER :: pPatchSrc, pPatchDes
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_region), POINTER :: pRegionSrc, pRegionDes
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_BufferDataSend.F90,v $ $Revision: 1.6 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_BufferDataSend',&
  'PLAG_BufferDataSend.F90' )
    
! get dimensions ------------------------------
  
  nCont = regions(iReg)%plagInput%nCont
  nCv   = CV_PLAG_LAST + nCont
  nDv   = DV_PLAG_LAST
  nTv   = TV_PLAG_LAST
  nAiv  = AIV_PLAG_LAST
  nArv  = ARV_PLAG_LAST
       
  nDimI = nAiv
  nDimR = 2*nArv +4*nCv 

! set pointer for Source Region ------------------------------
  
  pRegionSrc => regions(iReg)
      
  iLev     = pRegionSrc%currLevel
  nPatches = pRegionSrc%nPatches

  nDumCellsSrc = pRegionSrc%nDumCells

  pPlag => pRegionSrc%levels(iLev)%plag
  pCvPlagMass => pPlag%cvPlagMass

! loop over patches -----------------------------------------------------------

  DO iPatch = 1, nPatches

! - pointer is at Src region sending data to Des region -----------------------

    pPatchSrc => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = pPatchSrc%bcType
    iRegDes   = pPatchSrc%srcRegion
    iPatchDes = pPatchSrc%srcPatch

! - get patch information on source region ------------------------------------

    CALL RFLO_GetPatchIndices( pRegionSrc,pPatchSrc,iLev,ibegSrc,iendSrc, &
                               jbegSrc,jendSrc,kbegSrc,kendSrc )
    CALL RFLO_GetPatchDirection( pPatchSrc,idirSrc,jdirSrc,kdirSrc )
    CALL RFLO_GetCellOffset( pRegionSrc,iLev,iCOffSrc,ijCOffSrc )
                                  
! - region interface for various boundary conditions --------------------------

    IF ( (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
         (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
         (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
         (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
         (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) ) THEN

      IF ( regions(iRegDes)%procid /= global%myProcid ) THEN
        pRegionDes => regions(iRegDes)
        pPatchDes  => pRegionDes%levels(iLev)%patches(iPatchDes)

! -- set pointers and buffer size on source region ---------------------------
   
        pAivSrc    => pPatchSrc%bufferPlag%aiv
        pArvSrc    => pPatchSrc%bufferPlag%arv

        pCvSrc     => pPatchSrc%bufferPlag%cv
        pRhsSrc    => pPatchSrc%bufferPlag%rhs
        pRhsSumSrc => pPatchSrc%bufferPlag%rhsSum

        pArvOldSrc => pPatchSrc%bufferPlag%arvOld
        pCvOldSrc  => pPatchSrc%bufferPlag%cvOld
            
        pSendBuffI => pPatchSrc%bufferPlag%sendBuffI
        pSendBuffR => pPatchSrc%bufferPlag%sendBuffR        
        
        nBuffSizeSrc = pPatchSrc%bufferPlag%nBuffSize

! -- exit for null buffer size ------------------------------------------------
       
        IF ( nBuffSizeSrc == 0 ) GOTO 999
 
        nSendBuffI = nDimI * nBuffSizeSrc
        nSendBuffR = nDimR * nBuffSizeSrc
        
        pPatchSrc%bufferPlag%nSendBuffI = nSendBuffI
        pPatchSrc%bufferPlag%nSendBuffR = nSendBuffR
 
        iRequestPlag = pPatchSrc%bufferPlag%iRequest
        
        procDes = pRegionDes%procid

#ifdef PLAG_MPI_DEBUG                       
  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferSizeSend: iReg, nBuffSizeSrc, nDimI, nDimR, nSendBuffI, nSendBuffR, iRequestPlag  = ',&
                     iReg, nBuffSizeSrc, nDimI, nDimR, nSendBuffI, nSendBuffR, iRequestPlag   
#endif

! -- get patch information on destination region ------------------------------

        iLevDes = pRegionDes%currLevel
        nDumCellsDes = pRegionDes%nDumCells
        CALL RFLO_GetPatchIndices( pRegionDes,pPatchDes,iLevDes,ibegDes,iendDes, &
                                   jbegDes,jendDes,kbegDes,kendDes )
        CALL RFLO_GetPatchDirection( pPatchDes,idirDes,jdirDes,kdirDes )
        CALL RFLO_GetCellOffset( pRegionDes,iLevDes,iCOffDes,ijCOffDes ) 

#ifdef PLAG_MPI_DEBUG                       
  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferDataSend: iReg,ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc  = ',&
                     iReg,ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc

  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferDataSend: iRegDes,ibegDes,iendDes,jbegDes,jendDes,kbegDes,kendDes  = ',&
                     iRegDes,ibegDes,iendDes,jbegDes,jendDes,kbegDes,kendDes
#endif
                                        
! mapping between patches -----------------------------------------------------

                                                     l1DesDir =  1
        IF (pPatchSrc%srcL1beg > pPatchSrc%srcL1end) l1DesDir = -1
                                                     l2DesDir =  1
        IF (pPatchSrc%srcL2beg > pPatchSrc%srcL2end) l2DesDir = -1

        lbSrc    = pPatchSrc%lbound
        lbDes    = pPatchSrc%srcLbound
        alignSrc = pPatchSrc%align

        CALL RFLO_GetPatchMapping( lbSrc,lbDes,l1DesDir,l2DesDir,alignSrc, &
                                   idirSrc,jdirSrc,kdirSrc,                &
                                   idirDes,jdirDes,kdirDes,                &
                                   ibegSrc,iendSrc,jbegSrc,jendSrc,        &
                                   kbegSrc,kendSrc,                        &
                                   ibegDes,iendDes,jbegDes,jendDes,        &
                                   kbegDes,kendDes,mapMat )

! load send buffers while updating data ---------------------------------------  

!- Integer variables ----------------------------------------------------------
      
        DO iBuff = 1, nBuffSizeSrc 

! -- compute shift

          iBuffSend = nDimI*(iBuff-1) +1

!-- extract indices

          iCellsSrc = pAivSrc(AIV_PLAG_ICELLS,iBuff)          
          indexISrc = pAivSrc(AIV_PLAG_INDEXI,iBuff)
          indexJSrc = pAivSrc(AIV_PLAG_INDEXJ,iBuff)
          indexKSrc = pAivSrc(AIV_PLAG_INDEXK,iBuff)

!-- get mapping
   
          iCellsDes = IndIJKMap( indexISrc,indexJSrc,indexKSrc, &
                                 mapMat,iCOffDes,ijCOffDes )
          CALL GetIJK( iCellsDes,iCOffDes,ijCOffDes,nDumCellsDes,&
                       indexIDes,indexJDes,indexKDes)

!-- load integer data buffers

          iAiv = iBuffSend
          pSendBuffI(iAiv  ) = pAivSrc(AIV_PLAG_PIDINI,iBuff)
          pSendBuffI(iAiv+1) = pAivSrc(AIV_PLAG_REGINI,iBuff)
          pSendBuffI(iAiv+2) = iRegDes 
          pSendBuffI(iAiv+3) = iCellsDes
          pSendBuffI(iAiv+4) = indexIDes
          pSendBuffI(iAiv+5) = indexJDes
          pSendBuffI(iAiv+6) = indexKDes
          pSendBuffI(iAiv+7) = pAivSrc(AIV_PLAG_BURNSTAT,iBuff)
          pSendBuffI(iAiv+8) = pAivSrc(AIV_PLAG_STATUS,iBuff)

#ifdef PLAG_MPI_DEBUG                       
  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferDataSend-INT: procDes, iBuff, iAiv, iCellsDes, indexIDes, indexJDes, indexKDes = ',&
                     procDes, iBuff, iAiv, iCellsDes, indexIDes, indexJDes, indexKDes
#endif
 
        ENDDO ! iBuff

!- Real variables -------------------------------------------------------------
      
        DO iBuff = 1, nBuffSizeSrc 

! ---- compute shifts ---------------------------------------------------------

          iBuffSend = nDimR*(iBuff-1) +1
          iCv       = iBuffSend
          iRhs      = iBuffSend +nCv 
          iRhsSum   = iBuffSend +2*nCv 
          iCvOld    = iBuffSend +3*nCv  
          iArv      = iBuffSend +4*nCv 
          iArvOld   = iBuffSend +4*nCv +nArv

#ifdef PLAG_MPI_DEBUG     
  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferDataSend-REAL: procDes,iBuff, iBuffSend, iCv, iRhs, iRhsSum, iCvOld, iArv, iArvOld = ',&
                     procDes,iBuff, iBuffSend, iCv, iRhs, iRhsSum, iCvOld, iArv, iArvOld 
#endif
           
! ---- load real data buffers: cv ---------------------------------------------

          pSendBuffR(iCv  ) =  pCvSrc(CV_PLAG_XMOM,iBuff)
          pSendBuffR(iCv+1) =  pCvSrc(CV_PLAG_YMOM,iBuff) 
          pSendBuffR(iCv+2) =  pCvSrc(CV_PLAG_ZMOM,iBuff)
          pSendBuffR(iCv+3) =  pCvSrc(CV_PLAG_ENER,iBuff)
          pSendBuffR(iCv+4) =  pCvSrc(CV_PLAG_XPOS,iBuff)
          pSendBuffR(iCv+5) =  pCvSrc(CV_PLAG_YPOS,iBuff)
          pSendBuffR(iCv+6) =  pCvSrc(CV_PLAG_ZPOS,iBuff)
          pSendBuffR(iCv+7) =  pCvSrc(CV_PLAG_ENERVAPOR,iBuff)
          DO iCont = 1, nCont
            pSendBuffR(iCv+(CV_PLAG_LAST-1)+iCont) = pCvSrc(pCvPlagMass(iCont),iBuff)
          ENDDO ! iCont

! -- load real data buffers: rhs

          pSendBuffR(iRhs  ) =  pRhsSrc(CV_PLAG_XMOM,iBuff)
          pSendBuffR(iRhs+1) =  pRhsSrc(CV_PLAG_YMOM,iBuff) 
          pSendBuffR(iRhs+2) =  pRhsSrc(CV_PLAG_ZMOM,iBuff)
          pSendBuffR(iRhs+3) =  pRhsSrc(CV_PLAG_ENER,iBuff)
          pSendBuffR(iRhs+4) =  pRhsSrc(CV_PLAG_XPOS,iBuff)
          pSendBuffR(iRhs+5) =  pRhsSrc(CV_PLAG_YPOS,iBuff)
          pSendBuffR(iRhs+6) =  pRhsSrc(CV_PLAG_ZPOS,iBuff)
          pSendBuffR(iRhs+7) =  pRhsSrc(CV_PLAG_ENERVAPOR,iBuff)
          DO iCont = 1, nCont
            pSendBuffR(iRhs+(CV_PLAG_LAST-1)+iCont) = pRhsSrc(pCvPlagMass(iCont),iBuff)
          ENDDO ! iCont

!-- load real data buffers: rhsSum

          pSendBuffR(iRhsSum  ) =  pRhsSumSrc(CV_PLAG_XMOM,iBuff)
          pSendBuffR(iRhsSum+1) =  pRhsSumSrc(CV_PLAG_YMOM,iBuff) 
          pSendBuffR(iRhsSum+2) =  pRhsSumSrc(CV_PLAG_ZMOM,iBuff)
          pSendBuffR(iRhsSum+3) =  pRhsSumSrc(CV_PLAG_ENER,iBuff)
          pSendBuffR(iRhsSum+4) =  pRhsSumSrc(CV_PLAG_XPOS,iBuff)
          pSendBuffR(iRhsSum+5) =  pRhsSumSrc(CV_PLAG_YPOS,iBuff)
          pSendBuffR(iRhsSum+6) =  pRhsSumSrc(CV_PLAG_ZPOS,iBuff)
          pSendBuffR(iRhsSum+7) =  pRhsSumSrc(CV_PLAG_ENERVAPOR,iBuff)
          DO iCont = 1, nCont
            pSendBuffR(iRhsSum+(CV_PLAG_LAST-1)+iCont) = pRhsSumSrc(pCvPlagMass(iCont),iBuff)
          ENDDO ! iCont          

!-- load real data buffers: cvOld

          pSendBuffR(iCvOld  ) =  pCvOldSrc(CV_PLAG_XMOM,iBuff)
          pSendBuffR(iCvOld+1) =  pCvOldSrc(CV_PLAG_YMOM,iBuff) 
          pSendBuffR(iCvOld+2) =  pCvOldSrc(CV_PLAG_ZMOM,iBuff)
          pSendBuffR(iCvOld+3) =  pCvOldSrc(CV_PLAG_ENER,iBuff)
          pSendBuffR(iCvOld+4) =  pCvOldSrc(CV_PLAG_XPOS,iBuff)
          pSendBuffR(iCvOld+5) =  pCvOldSrc(CV_PLAG_YPOS,iBuff)
          pSendBuffR(iCvOld+6) =  pCvOldSrc(CV_PLAG_ZPOS,iBuff)
          pSendBuffR(iCvOld+7) =  pCvOldSrc(CV_PLAG_ENERVAPOR,iBuff)
          DO iCont = 1, nCont
            pSendBuffR(iCvOld+(CV_PLAG_LAST-1)+iCont) = pCvOldSrc(pCvPlagMass(iCont),iBuff)
          ENDDO ! iCont

!-- load real data buffers: arv  

          pSendBuffR(iArv  ) = pArvSrc(ARV_PLAG_SPLOAD,iBuff)
          pSendBuffR(iArv+1) = pArvSrc(ARV_PLAG_DISTOT,iBuff)

!-- load real data buffers: arvOld

          pSendBuffR(iArvOld  ) =  pArvOldSrc(ARV_PLAG_SPLOAD,iBuff)
          pSendBuffR(iArvOld+1) =  pArvOldSrc(ARV_PLAG_DISTOT,iBuff) 

        ENDDO ! iBuff 
       
#ifdef MPI

!--- integer
        
        tagDesI = regions(iRegDes)%localNumber &
                + PLAG_TAG_SHIFT +MPI_PATCHOFF*iPatchDes*iRegDes + procDes +1 

        IF(tagDesI .gt. global%mpiTagMax) tagDesI = MOD(tagDesI,global%mpiTagMax)
#ifdef PLAG_MPI_DEBUG                       
  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferDataSend-Integer: iReg, iRegDes, procDes, tagDesI = ',&
                     iReg, iRegDes, procDes,tagDesI
#endif
                                  
        CALL MPI_Isend( pSendBuffI,nSendBuffI,MPI_INTEGER,          &
                        procDes,tagDesI,global%mpiComm,             &
                        pPlag%requestsI(iRequestPlag),global%mpierr )
                        
        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ ) 

!--- real
      
        tagDesR = regions(iRegDes)%localNumber &
                + PLAG_TAG_SHIFT +MPI_PATCHOFF*iPatchDes*iRegDes + procDes +2 

        IF(tagDesR .gt. global%mpiTagMax) tagDesR = MOD(tagDesR,global%mpiTagMax)
#ifdef PLAG_MPI_DEBUG                       
  IF( nBuffSizeSrc /=0 )&
  WRITE(STDOUT,*) '  PLAG_BufferDataSend-Real: iReg, iRegDes, procDes, tagDesR = ',&
                     iReg, iRegDes, procDes,tagDesR
#endif
                                  
        CALL MPI_Isend( pSendBuffR,nSendBuffR,MPI_RFREAL,           &
                        procDes,tagDesR,global%mpiComm,             &
                        pPlag%requestsR(iRequestPlag),global%mpierr )
                        
        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ ) 

#endif

      ENDIF ! regions
    ENDIF ! bcType

999 CONTINUE

  ENDDO ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_BufferDataSend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_BufferDataSend.F90,v $
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
! Revision 1.2  2005/05/31 21:37:32  fnajjar
! Added ARV_PLAG_DISTOT for proper IO capabilities
!
! Revision 1.1  2004/12/01 20:56:56  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/04/09 23:04:12  fnajjar
! Added AIV_PLAG_STATUS to buffers being sent and received
!
! Revision 1.5  2004/03/21 00:43:32  fnajjar
! Fixed tags to be smaller number since Frost run-time system complains about size
!
! Revision 1.4  2004/03/12 23:42:31  fnajjar
! Bug fix for pSendBuffR with incorrect assignment
!
! Revision 1.3  2004/03/06 21:25:05  fnajjar
! Added PLAG_TAG_SHIFT to MPI-based communication tags
!
! Revision 1.2  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.1  2003/02/21 17:09:18  fnajjar
! Initial import
!
!******************************************************************************








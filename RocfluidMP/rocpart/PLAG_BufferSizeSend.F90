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
! Purpose: send buffer size to adjacent region on different processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: buffer size to other processors.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_BufferSizeSend.F90,v 1.4 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_BufferSizeSend( regions, iReg )

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
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef MPI
  INTEGER :: statusPlag(MPI_STATUS_SIZE)
#endif

  INTEGER :: bcType, iLev, iPatchDes,iRegDes, iRequestPlag,           &
             nArv, nAiv, nBuffI, nBuffR, nBuffSizeSrc, nCont,         &
             nCv, nDimBuffSize, nDv, nPatches,nSendBuffI,             &
             nSendBuffR, nTv, procDes, tagDes 
                                        
  TYPE(t_patch),  POINTER :: patchSrc, patchDes
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_BufferSizeSend.F90,v $ $Revision: 1.4 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_BufferSizeSend',&
  'PLAG_BufferSizeSend.F90' )
    
! get dimensions and set pointer ----------------------------------------------
  
  nCont = regions(iReg)%plagInput%nCont
  nCv   = CV_PLAG_LAST + nCont
  nDv   = DV_PLAG_LAST
  nTv   = TV_PLAG_LAST
  nAiv  = AIV_PLAG_LAST
  nArv  = ARV_PLAG_LAST
       
  nBuffI     = 2*nAiv
  nBuffR     = 2*nArv +4*nCv +nDv +nTv  
  nDimBuffSize = 1 
        
  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

  pPlag => regions(iReg)%levels(iLev)%plag
          
! loop over patches -----------------------------------------------------------

  DO iPatch = 1, nPatches

! - pointer is at Src region sending data to Des region -----------------------

    patchSrc => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patchSrc%bcType
    iRegDes   = patchSrc%srcRegion
    iPatchDes = patchSrc%srcPatch

! - region interface for various boundary conditions --------------------------

    IF ( (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
         (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
         (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
         (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
         (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) ) THEN

      IF ( regions(iRegDes)%procid /= global%myProcid ) THEN
        patchDes => regions(iRegDes)%levels(iLev)%patches(iPatchDes)
        
        nBuffSizeSrc = patchSrc%bufferPlag%nBuffSize
        
        nSendBuffI = nBuffI *nBuffSizeSrc
        nSendBuffR = nBuffR *nBuffSizeSrc
  
        patchSrc%bufferPlag%nSendBuffI = nSendBuffI
        patchSrc%bufferPlag%nSendBuffR = nSendBuffR 

        iRequestPlag = patchSrc%bufferPlag%iRequest

#ifdef PLAG_MPI_DEBUG              
  IF( nBuffSizeSrc /= 0 ) &
  WRITE(STDOUT,*) '  PLAG_BufferSizeSend: iReg, nBuffSizeSrc, nDimBuffSize, iRequestPlag  = ',&
                     iReg, nBuffSizeSrc, nDimBuffSize, iRequestPlag    
#endif

#ifdef MPI
        procDes = regions(iRegDes)%procid
        tagDes  = regions(iRegDes)%localNumber &
                + PLAG_TAG_SHIFT +MPI_PATCHOFF*iPatchDes*iRegDes + procDes 

        IF(tagDes .gt. global%mpiTagMax) tagDes = MOD(tagDes,global%mpiTagMax)

#ifdef PLAG_MPI_DEBUG          
  IF( nBuffSizeSrc /= 0 ) &
  WRITE(STDOUT,*) '  PLAG_BufferSizeSend: iReg, iRegDes, procDes, tagDes = ',&
                     iReg, iRegDes, procDes,tagDes
#endif
                                  
        CALL MPI_Isend( patchSrc%bufferPlag%nBuffSize,             &
                        nDimBuffSize,MPI_INTEGER,                  &
                        procDes,tagDes,global%mpiComm,             &
                        pPlag%requests(iRequestPlag),global%mpierr )
                        
        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
 
      ENDIF ! regions
    ENDIF ! bcType

  ENDDO ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_BufferSizeSend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_BufferSizeSend.F90,v $
! Revision 1.4  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:56:59  fnajjar
! Initial revision after changing case
!
! Revision 1.10  2004/04/09 23:05:54  fnajjar
! Added IF statement to activate only for non-null buffer size
!
! Revision 1.9  2004/03/21 00:43:32  fnajjar
! Fixed tags to be smaller number since Frost run-time system complains about size
!
! Revision 1.8  2004/03/06 21:25:05  fnajjar
! Added PLAG_TAG_SHIFT to MPI-based communication tags
!
! Revision 1.7  2003/05/27 19:14:16  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.6  2003/05/07 00:15:05  fnajjar
! Included I/O within ifdef PLAG_MPI_DEBUG construct
!
! Revision 1.5  2003/01/24 23:10:15  f-najjar
! Added to tagSrc and tagDes procId for Des region to avoid tag collision
!
! Revision 1.4  2003/01/24 22:37:55  f-najjar
! Made tagDes less prone to tag collision by multiplying by iRegDes
!
! Revision 1.3  2003/01/24 22:33:35  f-najjar
! Used generic ERR_NONE for MPI error trapping
!
! Revision 1.2  2003/01/23 17:51:34  f-najjar
! Add ModMPI for MPI communication
!
! Revision 1.1  2003/01/23 17:31:22  f-najjar
! Initial import for MPI
!
!******************************************************************************








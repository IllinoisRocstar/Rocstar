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
! Purpose: allocate memory for variables associated with buffer datastructure 
!          for all active regions on current processor.
!
! Description: none.
!
! Input: regions   = all regions
!        iReg      = region number.
!
! Output: region%levels%patch%buffPlag = Buffplag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_AllocateDataBuffers.F90,v 1.5 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_AllocateDataBuffers( regions, iReg )

  USE ModDataTypes 
  USE ModPartLag, ONLY    : t_plag_input, t_buffer_plag 
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_CECellsAllocateData
  USE ModMPI

#ifdef STATS
  USE PLAG_RFLO_ModStats, ONLY : PLAG_RFLO_CreateStatBuff
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
 
  INTEGER :: iReg

! ... loop variables
  INTEGER :: iLev, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, errorFlag, iRegDes, lbound, n1, n2, nAiv, nArv, &
             nBuffI, nBuffR,nBuffSizeI, nBuffSizeR, nBuffSizeTot,    &
             nCont, nCv, nDv, nPatchSize, nTv

  TYPE(t_patch),       POINTER :: pPatch
  TYPE(t_buffer_plag), POINTER :: pBuffPlag
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_AllocateDataBuffers.F90,v $ $Revision: 1.5 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global, 'PLAG_AllocateDataBuffers',&
  'PLAG_AllocateDataBuffers.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating Data Buffers for PLAG...'
  END IF ! global%verbLevel
  
! Get dimensions --------------------------------------------------------------

  nCont        = regions(iReg)%plagInput%nCont
  nBuffSizeTot = regions(iReg)%plagInput%nPclsBuffTot
          
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,regions(iReg)%nGridLevels

    pPlag => regions(iReg)%levels(iLev)%plag
    pPlag%nRequests = 0
   
    nAiv = pPlag%nAiv
    nArv = pPlag%nArv
            
    nCv  = pPlag%nCv
    nDv  = pPlag%nDv
    nTv  = pPlag%nTv  
    
    nBuffI = 2*nAiv
    nBuffR = 2*nArv +4*nCv +nDv +nTv
   
    DO iPatch=1,regions(iReg)%nPatches

      pPatch  => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType  =  pPatch%bcType
      lbound  =  pPatch%lbound
      iRegDes =  pPatch%srcRegion

      n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1
      n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
      nPatchSize  = n1*n2

      pBuffPlag => pPatch%bufferPlag

      IF ( (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
           (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
           (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) ) THEN 

        pBuffPlag%nBuffSizeTot = nBuffSizeTot

! - Allocate buffer data ------------------------------------------------------

        ALLOCATE( pBuffPlag%aiv(nAiv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%aiv' ) 
        END IF ! global%error 

        ALLOCATE( pBuffPlag%arv(nArv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%arv' )
        END IF ! global%error

        ALLOCATE( pBuffPlag%cv(nCv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%cv' ) 
        END IF ! global%error

        ALLOCATE( pBuffPlag%dv(nDv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%dv' ) 
        END IF ! global%error

        ALLOCATE( pBuffPlag%tv(nTv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%tv' ) 
        END IF ! global%error

        ALLOCATE( pBuffPlag%aivOld(nAiv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%aivOld' ) 
        END IF ! global%error 

        ALLOCATE( pBuffPlag%arvOld(nArv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%arvOld' ) 
        END IF ! global%error

        ALLOCATE( pBuffPlag%cvOld(nCv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%cvOld' ) 
        END IF ! global%error  

        ALLOCATE( pBuffPlag%rhs(nCv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%rhs' ) 
        END IF ! global%error  

        ALLOCATE( pBuffPlag%rhsSum(nCv,nBuffSizeTot),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%rhsSum' ) 
        END IF ! global%error 
          
! -- Initialize data --------------------------------------------------------
        
        pBuffPlag%nBuffSize    = 0
        pBuffPlag%nBuffSizeDes = 0

        pBuffPlag%aiv = 0
        pBuffPlag%arv = 0.0_RFREAL
        pBuffPlag%cv  = 0.0_RFREAL
        pBuffPlag%dv  = 0.0_RFREAL
        pBuffPlag%tv  = 0.0_RFREAL

        pBuffPlag%aivOld = 0.0_RFREAL
        pBuffPlag%arvOld = 0.0_RFREAL
        pBuffPlag%cvOld  = 0.0_RFREAL

        pBuffPlag%rhs    = 0.0_RFREAL
        pBuffPlag%rhsSum = 0.0_RFREAL

! - Allocate data for off-processor communication -----------------------------

        IF (regions(iRegDes)%procid /= global%myProcid) THEN  ! other processor 
          nBuffSizeI = nBuffI *nBuffSizeTot
          nBuffSizeR = nBuffR *nBuffSizeTot

          pBuffPlag%nSendBuffTotI = nBuffSizeI
          pBuffPlag%nSendBuffTotR = nBuffSizeR
 
          pBuffPlag%nRecvBuffTotI = nBuffSizeI
          pBuffPlag%nRecvBuffTotR = nBuffSizeR
 
          ALLOCATE( pBuffPlag%sendBuffR(nBuffSizeR),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%sendBuffR' ) 
          END IF ! global%error

          ALLOCATE( pBuffPlag%recvBuffR(nBuffSizeR),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%recvBuffR' ) 
          END IF ! global%error

          ALLOCATE( pBuffPlag%sendBuffI(nBuffSizeI),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%sendBuffI' ) 
          END IF ! global%error

          ALLOCATE( pBuffPlag%recvBuffI(nBuffSizeI),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%recvBuffI' ) 
          END IF ! global%error

! -- Initialize data --------------------------------------------------------

          pBuffPlag%sendBuffI = 0
          pBuffPlag%sendBuffR = 0.0_RFREAL
        
          pBuffPlag%recvBuffI = 0
          pBuffPlag%recvBuffR = 0.0_RFREAL

! -- Set MPI request  -------------------------------------------------------- 

          pPlag%nRequests = pPlag%nRequests + 1 
          pBuffPlag%iRequest  = pPlag%nRequests 
               
        ENDIF ! regions

      ELSE IF ((bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
               (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE)) THEN
        CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )  ! #### TEMPORARY ####

      ELSE 
        NULLIFY(pBuffPlag%aiv)
        NULLIFY(pBuffPlag%arv)
        NULLIFY(pBuffPlag%cv)
        NULLIFY(pBuffPlag%dv)
        NULLIFY(pBuffPlag%tv)
        NULLIFY(pBuffPlag%rhs)
        NULLIFY(pBuffPlag%rhsSum)
        NULLIFY(pBuffPlag%aivOld)
        NULLIFY(pBuffPlag%arvOld)
        NULLIFY(pBuffPlag%cvOld)
        NULLIFY(pBuffPlag%sendBuffR)
        NULLIFY(pBuffPlag%sendBuffI)
        NULLIFY(pBuffPlag%recvBuffR)
        NULLIFY(pBuffPlag%recvBuffI)
      ENDIF     ! bcType

    ENDDO             ! iPatch

! Allocate array for send requests --------------------------------------------
! Note: Need to take into account corners and edges

#ifdef MPI

!- request for buffer size

    ALLOCATE( pPlag%requests(pPlag%nRequests),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pPlag%requests' ) 
    END IF ! global%error

!- request for integer data buffer
    
    ALLOCATE( pPlag%requestsI(pPlag%nRequests),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pPlag%requestsI' ) 
    END IF ! global%error

!- request for real data buffer
 
    ALLOCATE( pPlag%requestsR(pPlag%nRequests),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pPlag%requestsR' ) 
    END IF ! global%error

#endif

  ENDDO   ! iLev

! Allocate pertinent data for corner and edge cells ---------------------------

  CALL PLAG_cECellsAllocateData(regions,iReg)

! Allocate buffer arrays for statistics ---------------------------------------

#ifdef STATS
  CALL PLAG_RFLO_CreateStatBuff(regions,iReg)
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_AllocateDataBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_AllocateDataBuffers.F90,v $
! Revision 1.5  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.2  2005/02/16 14:43:59  fnajjar
! Added call to allocate buffer arrays for statistics
!
! Revision 1.1  2004/12/01 20:56:49  fnajjar
! Initial revision after changing case
!
! Revision 1.12  2003/11/12 21:22:02  fnajjar
! Included Corner-Edge cells memory allocation
!
! Revision 1.11  2003/05/27 19:13:48  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.10  2003/02/06 16:14:46  f-najjar
! Included memory allocation for requests of integer and real data buffers
!
! Revision 1.9  2003/01/24 19:49:55  f-najjar
! Cleaned up definition of nAiv, nArv, nCv, nDv, nTv based on pPlag
!
! Revision 1.8  2003/01/24 19:40:23  f-najjar
! Bug fix to define nBuffL correctly
!
! Revision 1.7  2003/01/23 17:21:57  f-najjar
! Removed Hidden TABs
!
! Revision 1.6  2003/01/23 00:15:05  f-najjar
! Corrected buffer size for Real variables and defined nBuffI, nBuffL, nBuffR
!
! Revision 1.5  2003/01/23 00:01:02  f-najjar
! Redefined iRequest based on pBuffPlag
!
! Revision 1.4  2003/01/22 23:55:50  f-najjar
! Added MPI-related request data
!
! Revision 1.3  2003/01/16 22:35:14  f-najjar
! Activated buffers for on and off processor communication
!
! Revision 1.2  2003/01/13 19:04:52  f-najjar
! Added initialization for buffer data and renamed iRegSrc to iRegDes
!
! Revision 1.1  2002/10/25 14:13:59  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








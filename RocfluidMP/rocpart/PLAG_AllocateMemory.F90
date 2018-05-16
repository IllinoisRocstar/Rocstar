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
! Purpose: allocate memory for all variables associated with the Lagrangian 
!          particles (PLAG) for all active regions on current processor.
!
! Description: none.
!
! Input: region = current region
!
! Output: region%plag = plag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_AllocateMemory.F90,v 1.5 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_AllocateMemory( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetNodeOffset
  USE PLAG_ModInterfaces, ONLY : PLAG_AllocateMemoryTile, &
                                 PLAG_InitMemory,         &
                                 PLAG_InjcTileInitialize
  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters

#ifdef STATS
  USE PLAG_ModEulerian, ONLY: PLAG_CreateEulerianField
  USE PLAG_ModStats, ONLY: PLAG_CreateStat
#endif

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iCont, iLev

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: errorFlag, ibn, idnbeg, idnend, ien, iNOff, ijNOff, &
             jdnbeg, jdnend, kdnbeg, kdnend, nAiv, nArv, nCont,  &
             nCv, nDv, nEv, nPclsMax, nTv, maxDisEdges

  TYPE(t_level),        POINTER :: pLevel
  TYPE(t_plag),         POINTER :: pPlag
  TYPE(t_global),       POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_AllocateMemory.F90,v $ $Revision: 1.5 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_AllocateMemory',&
  'PLAG_AllocateMemory.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating memory for PLAG...'
  END IF ! global%verbLevel

! Set max discrete edges ------------------------------------------------------
  
  maxDisEdges = region%inrtInput%maxDisEdges 

! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    pLevel    => region%levels(iLev)
    pPlag     => region%levels(iLev)%plag

! - Get particle dimensions ---------------------------------------------------

    nPclsMax = region%plagInput%nPclsMax
    nCont    = region%plagInput%nCont

! - Set pointers --------------------------------------------------------------

    nAiv = pPlag%nAiv
    nArv = pPlag%nArv
            
    nCv  = pPlag%nCv
    nDv  = pPlag%nDv
    nTv  = pPlag%nTv
    nEv  = pPlag%nEv

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,1015) '   nPclsMax',nPclsMax
      WRITE(STDOUT,1010) '   nCont   ',nCont
      WRITE(STDOUT,1010) '   nAiv    ',nAiv
      WRITE(STDOUT,1010) '   nArv    ',nArv
      WRITE(STDOUT,1010) '   nCv     ',nCv
      WRITE(STDOUT,1010) '   nDv     ',nDv
      WRITE(STDOUT,1010) '   nTv     ',nTv
      WRITE(STDOUT,1010) '   nEv     ',nEv
    END IF ! global%verbLevel
     
! - Lagrangian particles variables --------------------------------------------

    ALLOCATE( pPlag%aiv(nAiv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%aiv' ) 
    END IF ! global%error

    ALLOCATE( pPlag%arv(nArv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%arv' ) 
    END IF ! global%error
    
    ALLOCATE( pPlag%cv (nCv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%cv' ) 
    END IF ! global%error

    ALLOCATE( pPlag%dv (nDv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%dv' ) 
    END IF ! global%error
    
    ALLOCATE( pPlag%tv (nTv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%tv' ) 
    END IF ! global%error
    
    ALLOCATE( pPlag%rhs   (nCv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%rhs' ) 
    END IF ! global%error
          
    ALLOCATE( pPlag%rhsSum(nCv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%rhsSum' ) 
    END IF ! global%error

    ALLOCATE( pPlag%aivOld(nAiv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%aivOld' ) 
    END IF ! global%error

    ALLOCATE( pPlag%arvOld(nArv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%arvOld' ) 
    END IF ! global%error
    
    ALLOCATE( pPlag%cvOld (nCv,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%cvOld' ) 
    END IF ! global%error

! - Lagrangian particles mass and volume indices ------------------------------
          
    ALLOCATE( pPlag%cvPlagMass(nCont),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%cvPlagMass' )
    END IF ! global%error
          
    ALLOCATE( pPlag%dvPlagVolu(nCont),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%dvPlagVolu' )
    END IF ! global%error
 
    DO iCont = 1, nCont
      pPlag%cvPlagMass(iCont) = CV_PLAG_LAST  +iCont
      pPlag%dvPlagVolu(iCont) = DV_PLAG_LAST  +iCont
    END DO ! iCont

! - Get dimensions and pointers -----------------------------------------------

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

! - Allocate array for face centroids ----------------------------------------- 
   
    ALLOCATE( pPlag%fc(ZCOORD,KCOORD,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%fc' )
    END IF ! global%error

! - Allocate array for face vectors ------------------------------------------- 

    ALLOCATE( pPlag%si(ZCOORD,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%si' )
    END IF ! global%error

    ALLOCATE( pPlag%sj(ZCOORD,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%sj' )
    END IF ! global%error

    ALLOCATE( pPlag%sk(ZCOORD,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'pPlag%sk' )
    END IF ! global%error

! - Allocate source for INRT ----------------------------------------------------

!    CALL PLAG_INRT_AllocMem( region, pPlag ) 

    IF ( maxDisEdges > 0 ) THEN
      ALLOCATE( pPlag%inrtSources(maxDisEdges,nPclsMax), stat=errorFlag )

      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'pPlag%inrtSources'  )
      END IF ! global%error

    ELSE
      NULLIFY( pPlag%inrtSources )

    END IF ! maxDisEdges

! - Allocate arrays for statistics --------------------------------------------

#ifdef STATS
    IF ( ( global%flowType == FLOW_UNSTEADY ) .AND. &
         ( global%doStat == ACTIVE )                ) THEN
      CALL PLAG_CreateEulerianField( region, pPlag )
      CALL PLAG_CreateStat( region, pPlag )
    END IF ! global%flowType
#endif

  ENDDO   ! iLev

! Allocate memory for tiles ---------------------------------------------------

  CALL PLAG_AllocateMemoryTile( region )

! Initialize memory -----------------------------------------------------------

  CALL PLAG_InitMemory( region )

! Initialize tile datastructure -----------------------------------------------

  CALL PLAG_InjcTileInitialize( region ) 

! finalize

  CALL DeregisterFunction( global )
  
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I7)

END SUBROUTINE PLAG_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_AllocateMemory.F90,v $
! Revision 1.5  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/06 23:13:13  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.2  2005/01/08 20:41:38  fnajjar
! Included memory allocation for PLAG statistics
!
! Revision 1.1  2004/12/01 20:56:50  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/07/26 17:05:51  fnajjar
! moved allocation of inrtSources into Rocpart
!
! Revision 1.6  2004/03/01 21:57:54  fnajjar
! Removed definitions of nAiv, nCv, nDv, nTv since now set in PLAG_derivedInputValues
!
! Revision 1.5  2003/11/03 21:19:29  fnajjar
! Added face vectors to PLAG infrastructure
!
! Revision 1.4  2003/05/27 19:14:16  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.3  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.2  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
! Revision 1.1  2002/10/25 14:13:59  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








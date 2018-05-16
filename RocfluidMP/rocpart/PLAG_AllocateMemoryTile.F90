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
! Purpose: allocate memory for variables associated with Tile datastructure 
!          for current region.
!
! Description: none.
!
! Input: region = current region
!
! Output: region%levels%patch%tilePlag = Tileplag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_AllocateMemoryTile.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_AllocateMemoryTile( region )

  USE ModDataTypes 
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag 
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iCont, iLev, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, errorFlag, n1, n2, nCont, nCv, nDv, nTile

  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global),    POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_AllocateMemoryTile.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global, 'PLAG_AllocateMemoryTile',&
  'PLAG_AllocateMemoryTile.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating tile memory for PLAG...'
  END IF ! global%verbLevel
  
! loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels
    DO iPatch=1,region%nPatches
  
      pPatch => region%levels(iLev)%patches(iPatch)
      bcType =  pPatch%bcType
      nCont  =  region%plagInput%nCont

      IF ( bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE ) THEN
         n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1    
         n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
         nTile   = n1*n2

         pTilePlag => pPatch%tilePlag
      
         pTilePlag%nCv = CV_TILE_LAST + nCont
         pTilePlag%nDv = DV_TILE_LAST

         nCv  = pTilePlag%nCv
         nDv  = pTilePlag%nDv

         ALLOCATE( pTilePlag%nPclsInjc(nTile),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%nPclsInjc' ) 
         END IF ! global%error
         
         ALLOCATE( pTilePlag%cv(nCv,nTile),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%cv' ) 
         END IF ! global%error

         ALLOCATE( pTilePlag%cvOld(nCv,nTile),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%cvOld' ) 
         END IF ! global%error
          
         ALLOCATE( pTilePlag%dv(nDv,nTile),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%dv' ) 
         END IF ! global%error

         ALLOCATE( pTilePlag%rhs(nCv,nTile),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%rhs' ) 
         END IF ! global%error

         ALLOCATE( pTilePlag%rhsSum(nCv,nTile),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%rhsSum' ) 
         END IF ! global%error

         ALLOCATE( pTilePlag%cvTileMass(nCont),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pTilePlag%cvTileMass' ) 
         END IF ! global%error

         DO iCont = 1, nCont
           pTilePlag%cvTileMass(iCont) = CV_TILE_LAST +iCont
         END DO ! iCont

      ENDIF  ! bcType

    ENDDO    ! iPatch

  ENDDO   ! iLev

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_AllocateMemoryTile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_AllocateMemoryTile.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:56:51  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2002/10/25 14:13:59  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








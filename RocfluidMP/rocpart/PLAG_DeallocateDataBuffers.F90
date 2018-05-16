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
! Purpose: deallocate memory for variables associated with buffer datastructure 
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
! $Id: PLAG_DeallocateDataBuffers.F90,v 1.4 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_DeallocateDataBuffers( regions, iReg )

  USE ModDataTypes 
  USE ModPartLag, ONLY    : t_buffer_plag 
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
 
  INTEGER :: iReg
   
! ... loop variables
  INTEGER :: iLev, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, errorFlag, iRegSrc, lbound, n1, n2, nCont, &
             nPatchSize, nBuffSizeTot

  TYPE(t_patch),       POINTER :: pPatch
  TYPE(t_buffer_plag), POINTER :: pBuffPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_DeallocateDataBuffers.F90,v $ $Revision: 1.4 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global, 'PLAG_DeallocateDataBuffers',&
  'PLAG_DeallocateDataBuffers.F90' )

! Get dimensions --------------------------------------------------------------

  nBuffSizeTot = regions(iReg)%plagInput%nPclsBuffTot

! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,regions(iReg)%nGridLevels
    DO iPatch=1,regions(iReg)%nPatches
  
      pPatch  => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType  =  pPatch%bcType
      lbound  =  pPatch%lbound
      iRegSrc =  pPatch%srcRegion
      
      n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1    
      n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
      nPatchSize  = n1*n2

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN 

        pBuffPlag => pPatch%bufferPlag
      
        pBuffPlag%nBuffSizeTot = nBuffSizeTot

! - Deallocate data for same-processor communication -------------------------- 
      
        IF (regions(iRegSrc)%procid == global%myProcid) THEN 
          DEALLOCATE( pBuffPlag%aiv,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%aiv' ) 
          END IF ! global%error
          
          DEALLOCATE( pBuffPlag%arv,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%arv' ) 
          END IF ! global%error

          DEALLOCATE( pBuffPlag%cv,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%cv' ) 
          END IF ! global%error

          DEALLOCATE( pBuffPlag%dv,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%dv' ) 
          END IF ! global%error

          DEALLOCATE( pBuffPlag%tv,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%tv' ) 
          END IF ! global%error
          
          DEALLOCATE( pBuffPlag%aivOld,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%aivOld' ) 
          END IF ! global%error
          
          DEALLOCATE( pBuffPlag%arvOld,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%arvOld' ) 
          END IF ! global%error

          DEALLOCATE( pBuffPlag%cvOld,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%cvOld' ) 
          END IF ! global%error
          DEALLOCATE( pBuffPlag%rhs,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%rhs' ) 
          END IF ! global%error

          DEALLOCATE( pBuffPlag%rhsSum,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%rhsSum' ) 
          END IF ! global%error
          
        END IF !regions

! - Allocate data for off-processor communication ----------------------------- 

        IF (regions(iRegSrc)%procid /= global%myProcid) THEN  ! other processor        

          DEALLOCATE( pBuffPlag%sendBuffR,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%sendBuffR' ) 
          END IF ! global%error
     
          DEALLOCATE( pBuffPlag%recvBuffR,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%recvBuffR' )
          END IF ! global%error

          DEALLOCATE( pBuffPlag%sendBuffI,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%sendBuffI' ) 
          END IF ! global%error 
     
          DEALLOCATE( pBuffPlag%recvBuffI,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pBuffPlag%recvBuffI' ) 
          END IF ! global%error

        ENDIF ! regions

      ELSE IF ((bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
               (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE)) THEN
        CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )  ! #### TEMPORARY ####
      ENDIF     ! bcType
                     
    ENDDO             ! iPatch

  ENDDO   ! iLev

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_DeallocateDataBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_DeallocateDataBuffers.F90,v $
! Revision 1.4  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:29  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2003/05/27 19:19:57  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.1  2002/10/25 14:15:43  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








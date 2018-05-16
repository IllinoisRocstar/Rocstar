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
! Purpose: deallocate memory for main variables associated with the Lagrangian 
!          particles (PLAG) for all active regions on current processor.
!
! Description: none.
!
! Input: region = current region data
!        iReg   = current region
!
! Output: region%plag = plag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_DeallocateMemoryPost.F90,v 1.3 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_DeallocateMemoryPost( region, iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

  INTEGER, INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: iCont, iLev

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: errorFlag

  TYPE(t_level),        POINTER :: pLevel
  TYPE(t_plag),         POINTER :: pPlag
  TYPE(t_global),       POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_DeallocateMemoryPost.F90,v $ $Revision: 1.3 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_DeallocateMemoryPost',&
  'PLAG_DeallocateMemoryPost.F90' )

  IF ( iReg == 1 .AND. global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Deallocating memory for PLAG...'
  END IF ! global%verbLevel
  
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    pLevel    => region%levels(iLev)
    pPlag     => region%levels(iLev)%plag

! - Lagrangian particles variables --------------------------------------------

    DEALLOCATE( pPlag%aiv,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%aiv' ) 
    END IF ! global%error

    DEALLOCATE( pPlag%arv,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%arv' ) 
    END IF ! global%error
    
    DEALLOCATE( pPlag%cv,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%cv' ) 
    END IF ! global%error

    DEALLOCATE( pPlag%dv,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%dv' ) 
    END IF ! global%error
    
    DEALLOCATE( pPlag%tv,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%tv' ) 
    END IF ! global%error
    
! - Lagrangian particles mass and volume indices ------------------------------
          
    DEALLOCATE( pPlag%cvPlagMass,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%cvPlagMass' )
    END IF ! global%error
          
    DEALLOCATE( pPlag%dvPlagVolu,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'pPlag%dvPlagVolu' )
    END IF ! global%error

  ENDDO   ! iLev

! finalize

  IF ( iReg == global%nRegions .AND. global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Deallocating memory for PLAG done...'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_DeallocateMemoryPost

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_DeallocateMemoryPost.F90,v $
! Revision 1.3  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/02/16 14:52:40  fnajjar
! Initial import
!
!******************************************************************************








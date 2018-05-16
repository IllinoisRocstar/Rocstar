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
! Purpose: allocate memory for main variables associated with the Lagrangian 
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
! $Id: PLAG_AllocateMemoryPost.F90,v 1.5 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_AllocateMemoryPost( region, iReg )

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

  INTEGER :: errorFlag, kdnend, nAiv, nArv, nCont, nCv, nDv, nPclsMax, nTv

  TYPE(t_level),        POINTER :: pLevel
  TYPE(t_plag),         POINTER :: pPlag
  TYPE(t_global),       POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_AllocateMemoryPost.F90,v $ $Revision: 1.5 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_AllocateMemoryPost',&
  'PLAG_AllocateMemoryPost.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating memory for PLAG...'
  END IF ! global%verbLevel
  
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    pLevel    => region%levels(iLev)
    pPlag     => region%levels(iLev)%plag

! - Get particle dimensions ---------------------------------------------------

    nPclsMax = region%plagInput%nPclsMax
    nCont    = region%plagInput%nCont

    IF ( iReg == 1 ) &
      PRINT*,' PLAG_AllocateMemoryPost: nPclsMax = ',nPclsMax

! - Set pointers --------------------------------------------------------------
        
    pPlag%nCv  = CV_PLAG_LAST + nCont
    pPlag%nDv  = DV_PLAG_LAST + nCont
    pPlag%nTv  = TV_PLAG_LAST
    pPlag%nAiv = AIV_PLAG_LAST
    pPlag%nArv = ARV_PLAG_LAST
    
    nAiv = pPlag%nAiv
    nArv = pPlag%nArv
            
    nCv  = pPlag%nCv
    nDv  = pPlag%nDv
    nTv  = pPlag%nTv

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,1015) '   nPclsMax',nPclsMax
      WRITE(STDOUT,1010) '   nCont   ',nCont
      WRITE(STDOUT,1010) '   nAiv    ',nAiv
      WRITE(STDOUT,1010) '   nArv    ',nArv
      WRITE(STDOUT,1010) '   nCv     ',nCv
      WRITE(STDOUT,1010) '   nDv     ',nDv
      WRITE(STDOUT,1010) '   nTv     ',nTv
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

  ENDDO   ! iLev

! finalize

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating memory for PLAG done...'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )
  
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I7)

END SUBROUTINE PLAG_AllocateMemoryPost

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_AllocateMemoryPost.F90,v $
! Revision 1.5  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/06 23:27:44  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.2  2005/02/16 14:50:19  fnajjar
! Cosmetic change to include verbosity at end of routine
!
! Revision 1.1  2004/12/01 22:00:46  fnajjar
! Initial revision after changing case
!
! Revision 1.1.1.1  2003/05/06 16:14:38  fnajjar
! Import of postprocessing tool for Rocpart
!
!******************************************************************************








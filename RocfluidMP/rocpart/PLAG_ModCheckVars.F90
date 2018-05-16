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
! ******************************************************************************
!
! Purpose: Suite of routines to check validity and positivy of PLAG variables.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModCheckVars.F90,v 1.7 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModCheckVars

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag
  USE ModMPI
  
  USE PLAG_ModParameters
  
  USE ModTools, ONLY: IsNan

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_CheckPositivity, &
            PLAG_CheckValidity
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: PLAG_ModCheckVars.F90,v $ $Revision: 1.7 $' 
                      
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
! Purpose: Check for posivity of variables
!
! Description: None.
!
! Input: 
!   pRegion     Region data
!
! Output: None.
!
! Notes: 
!  1. dv field has not been updated when this check is done. However,
!     the trap will be quite helpful.
!
! ******************************************************************************

  SUBROUTINE PLAG_CheckPositivity(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString

    INTEGER, PARAMETER :: MAX_NEGATIVE_LOCS = 10
    INTEGER :: icg,idIni,iPcl,nLocs,nPcls,regIni,stat
#ifdef RFLO
    INTEGER :: iLev
#endif
    INTEGER :: loc(MAX_NEGATIVE_LOCS,MIN_VAL:MAX_VAL)
    INTEGER, DIMENSION(:,:), POINTER :: pAiv

    REAL(RFREAL) :: diam,ener,mass,temp,xpos,ypos,zpos
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv
  
    TYPE(t_global), POINTER :: global
    TYPE(t_plag),   POINTER :: pPlag

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_CheckPositivity',&
  'PLAG_ModCheckVars.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

#ifdef RFLO
    iLev  =  pRegion%currLevel
    pPlag => pRegion%levels(iLev)%plag
#endif

#ifdef RFLU
    pPlag => pRegion%plag
#endif 

    pCv => pPlag%cv
    pDv => pPlag%dv
    pAiv => pPlag%aiv

    nPcls = pPlag%nPcls 
    nLocs = 0

! ******************************************************************************
!   Loop over particles and check for positivity
! ******************************************************************************

    DO iPcl = 1,nPcls
      mass  = SUM( pCv(pPlag%cvPlagMass(:),iPcl) )
      ener  = pCv(CV_PLAG_ENER,iPcl)
      diam  = pDv(DV_PLAG_DIAM,iPcl)
      temp  = pDv(DV_PLAG_TEMP,iPcl)
      xpos  = pCv(CV_PLAG_XPOS,iPcl)
      ypos  = pCv(CV_PLAG_YPOS,iPcl)
      zpos  = pCv(CV_PLAG_ZPOS,iPcl)
      
      icg    = pAiv(AIV_PLAG_ICELLS,iPcl)
      idIni  = pAiv(AIV_PLAG_PIDINI,iPcl)
      regIni = pAiv(AIV_PLAG_REGINI,iPcl)
      stat   = pAiv(AIV_PLAG_STATUS,iPcl)

      IF ( stat /= PLAG_STATUS_KEEP) CYCLE

      IF ( (mass <= 0.0_RFREAL) .OR. (ener <= 0.0_RFREAL) ) THEN 
        nLocs = nLocs + 1   

        IF ( nLocs == 1 ) THEN 
          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                  'Negative positive-definite variables detected!'
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Lagrangian Particle (PLAG).'        
          
          WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                              global%currentTime              
                                            
          WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                           pRegion%iRegionGlobal  
          WRITE(STDOUT,'(A,6X,A,11(1X,A))') SOLVER_NAME,'#', &
                                           '    iPcl     ', &
                                           '    idIni    ', &
                                           '    RegIni   ', &
                                           '    icg      ', &
                                           '    Mass     ', &
                                           '  x-location ', &
                                           '  y-location ', &
                                           '  z-Location ', &
                                           '   Energy    ', &
                                           '   Diameter  '       
        END IF ! nLocs

        IF ( nLocs <= MAX_NEGATIVE_LOCS ) THEN 
          WRITE(STDOUT,'(A,4X,4(1X,I8),6(1X,E13.6))') SOLVER_NAME,iPcl,    & 
                                                      idIni,regIni,icg,    &
                                                      mass,xpos,ypos,zpos, &
                                                      ener,diam                                
          loc(nLocs,MIN_VAL:MAX_VAL) = iPcl
        END IF ! nLocs
      END IF ! cv    
    END DO ! iPcl

! ******************************************************************************
!   Write out message and call error handling routine
! ******************************************************************************

    IF ( nLocs > 0 ) THEN 
      IF ( nLocs > MAX_NEGATIVE_LOCS ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
              'Only wrote the first',MAX_NEGATIVE_LOCS,'of',nLocs, & 
              'particles with negative positive-definite variables.'
      END IF ! nLocs
    
      CALL ErrorStop(global,ERR_NEGATIVE_POSDEF,__LINE__)   
    END IF ! nLocs

! ******************************************************************************
!   End
! ******************************************************************************

   CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CheckPositivity






! ******************************************************************************
! Purpose: Check validity of variables
!
! Description: None.
!
! Input: 
!   pRegion     Region data
!
! Output: None.
!
! Notes: 
!  1. dv field has not been updated when this check is done. However,
!     the trap will be quite helpful.
!
! ******************************************************************************

  SUBROUTINE PLAG_CheckValidity(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString
    INTEGER, PARAMETER :: MAX_INVALID_LOCS = 10
    INTEGER :: icg,idIni,iPcl,nLocs,nPcls,regIni,stat
#ifdef RFLO
    INTEGER :: iLev
#endif
    INTEGER :: loc(MAX_INVALID_LOCS,MIN_VAL:MAX_VAL)
    INTEGER, DIMENSION(:,:), POINTER :: pAiv

    REAL(RFREAL) :: diam,ener,mass,temp,xmom,xpos,ymom,ypos,zmom,zpos
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv
  
    TYPE(t_global), POINTER :: global
    TYPE(t_plag),   POINTER :: pPlag

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_CheckValidity',&
  'PLAG_ModCheckVars.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

#ifdef RFLO                                
    iLev  =  pRegion%currLevel               
    pPlag => pRegion%levels(iLev)%plag       
#endif                                     
                                           
#ifdef RFLU                                
    pPlag => pRegion%plag                  
#endif                                     
        
    pCv => pPlag%cv
    pDv => pPlag%dv
    pAiv => pPlag%aiv
                                                      
    nPcls = pPlag%nPcls 
    nLocs = 0

! ******************************************************************************
!   Loop over particles and check for validity
! ******************************************************************************

    DO iPcl = 1,nPcls
      mass  = SUM( pCv(pPlag%cvPlagMass(:),iPcl) )
      xmom  = pCv(CV_PLAG_XMOM,iPcl)
      ymom  = pCv(CV_PLAG_YMOM,iPcl)
      zmom  = pCv(CV_PLAG_ZMOM,iPcl)
      ener  = pCv(CV_PLAG_ENER,iPcl)
      diam  = pDv(DV_PLAG_DIAM,iPcl)
      temp  = pDv(DV_PLAG_TEMP,iPcl)
      xpos  = pCv(CV_PLAG_XPOS,iPcl)
      ypos  = pCv(CV_PLAG_YPOS,iPcl)
      zpos  = pCv(CV_PLAG_ZPOS,iPcl)
      
      icg    = pAiv(AIV_PLAG_ICELLS,iPcl)
      idIni  = pAiv(AIV_PLAG_PIDINI,iPcl)
      regIni = pAiv(AIV_PLAG_REGINI,iPcl)
      stat   = pAiv(AIV_PLAG_STATUS,iPcl)

      IF ( stat /= PLAG_STATUS_KEEP) CYCLE

      IF ( (IsNan(mass) .EQV. .TRUE.) .OR. & 
           (IsNan(xmom) .EQV. .TRUE.) .OR. &
           (IsNan(ymom) .EQV. .TRUE.) .OR. & 
           (IsNan(zmom) .EQV. .TRUE.) .OR. & 
           (IsNan(ener) .EQV. .TRUE.) .OR. &
           (IsNan(xpos) .EQV. .TRUE.) .OR. &
           (IsNan(ypos) .EQV. .TRUE.) .OR. &
           (IsNan(zpos) .EQV. .TRUE.) .OR. &
           (IsNan(diam) .EQV. .TRUE.) .OR. &
           (IsNan(temp) .EQV. .TRUE.) ) THEN
        nLocs = nLocs + 1   

        IF ( nLocs == 1 ) THEN 
          WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, & 
                'Invalid variables detected!'
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Lagrangian Particle (PLAG).' 
          WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                              global%currentTime              

          WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                           pRegion%iRegionGlobal  
          WRITE(STDOUT,'(A,6X,A,11(1X,A))') SOLVER_NAME,'#', &
                                           '    iPcl     ', &
                                           '    idIni    ', &
                                           '    RegIni   ', &
                                           '    icg      ', &
                                           '    Mass     ', &
                                           '  x-location ', &
                                           '  y-location ', &
                                           '  z-Location ', &
                                           '   Energy    ', &
                                           '   Diameter  '       
        END IF ! nLocs

        IF ( nLocs <= MAX_INVALID_LOCS ) THEN 
          WRITE(STDOUT,'(A,4X,4(1X,I8),6(1X,E13.6))') SOLVER_NAME,iPcl,    & 
                                                      idIni,regIni,icg,    &
                                                      mass,xpos,ypos,zpos, &
                                                      ener,diam                                
          loc(nLocs,MIN_VAL:MAX_VAL) = iPcl
        END IF ! nLocs
      END IF ! cv    
    END DO ! iPcl

! ******************************************************************************
!   Write out message and call error handling routine
! ******************************************************************************

    IF ( nLocs > 0 ) THEN 
      IF ( nLocs > MAX_INVALID_LOCS ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
              'Only wrote the first',MAX_INVALID_LOCS,'of',nLocs, & 
              'particles with invalid variables.'    
      END IF ! nLocs
    
      CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)   
    END IF ! nLocs

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_CheckValidity








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModCheckVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModCheckVars.F90,v $
! Revision 1.7  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.4  2005/12/19 16:49:28  fnajjar
! Added if statement to check validity for kept particles else cycle
!
! Revision 1.3  2005/12/07 20:05:25  fnajjar
! Removed check on diam and temp as vals not updated yet when using injection bc
!
! Revision 1.2  2005/12/05 19:28:38  fnajjar
! Bug fix for PLAG pointers in RFLO
!
! Revision 1.1  2005/12/01 21:53:48  fnajjar
! Initial version
!
! ******************************************************************************









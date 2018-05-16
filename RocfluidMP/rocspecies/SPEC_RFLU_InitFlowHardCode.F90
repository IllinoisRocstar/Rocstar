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
! Purpose: Initialize species flow field in a region using hard-coded values.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Initialize state vector in primitive form. This is done so can compute
!      gas properties based on species mass fractions.
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_InitFlowHardCode.F90,v 1.17 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_InitFlowHardCode(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
             
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,RCSIdentString
  INTEGER :: icg,iSpec
  REAL(RFREAL) :: x,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: SPEC_RFLU_InitFlowHardCode.F90,v $ $Revision: 1.17 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_InitFlowHardCode', &
                        'SPEC_RFLU_InitFlowHardCode.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing species flow field from hard code...'

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid      => pRegion%grid
  pCvMixt    => pRegion%mixt%cv
  pCvSpec    => pRegion%spec%cv
  pMixtInput => pRegion%mixtInput
  
  pRegion%spec%cvState = CV_MIXT_STATE_PRIM   
  
! ******************************************************************************
! Initialize flow field based on user input
! ******************************************************************************

  SELECT CASE ( global%casename )

! ==============================================================================
!   Generic compressible gravity current
! ==============================================================================
  
    CASE ( "gcgc" )
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies

      DO icg = 1,pGrid%nCellsTot           
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)

        IF ( (x < pMixtInput%prepRealVal1) .AND. & 
             (y < pMixtInput%prepRealVal2) ) THEN 
          pCvSpec(1,icg) = 0.04_RFREAL
          pCvSpec(2,icg) = 0.96_RFREAL                        
        ELSE 
          pCvSpec(1,icg) = 1.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL            
        END IF ! x
      END DO ! icg                       

! ==============================================================================
!   Generic multiphase jet
! ==============================================================================
  
    CASE ( "gmpjet" )
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies
    
      DO icg = 1,pGrid%nCellsTot           
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < 0.0_RFREAL ) THEN 
          pCvSpec(1,icg) = 0.5_RFREAL
          pCvSpec(2,icg) = 0.5_RFREAL                        
        ELSE 
          pCvSpec(1,icg) = 1.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL            
        END IF ! x
      END DO ! icg                       

! ==============================================================================
!   Kieffer jet 
! ==============================================================================

    CASE ( "kjet2v3mp","kjet2v4mp","kjet2v5mp" )
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies

      DO icg = 1,pGrid%nCellsTot     
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < pMixtInput%prepRealVal8 ) THEN 
          pCvSpec(1,icg) = pMixtInput%prepRealVal9
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal9
        ELSE 
          pCvSpec(1,icg) = pMixtInput%prepRealVal10
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal10
        END IF ! x
      END DO ! icg      

! ==============================================================================
!   Multiphase Shocktube: Water-Air
! ==============================================================================

    CASE ( "MShock_H2O_Air001", &
           "MShock_H2O_Air002"  )
      IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization only valid with gas-liq model.')
      END IF ! pRegion%mixtInput%gasModel   

      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies    

      DO icg = 1,pGrid%nCellsTot
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < 0.7_RFREAL ) THEN
          pCvSpec(1,icg) = 0.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE
          pCvSpec(1,icg) = 1.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        END IF ! x
      END DO ! icg   

! ==============================================================================
!   Multiphase Shocktube: Air-Air-He
! ==============================================================================
      
    CASE ( "MShock_Air_Air_He001" )
      IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization only valid with gas-liq model.')
      END IF ! pRegion%mixtInput%gasModel   

      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies    

      DO icg = 1,pGrid%nCellsTot
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < 0.25_RFREAL ) THEN
          pCvSpec(1,icg) = 1.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE IF ( x > 0.25_RFREAL .AND. x < 0.5_RFREAL ) THEN
          pCvSpec(1,icg) = 1.0_RFREAL 
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE
          pCvSpec(1,icg) = 0.0_RFREAL 
          pCvSpec(2,icg) = 1.0_RFREAL 
        END IF ! x
      END DO ! icg


! ============================================================================== 
!   Nozzle cavity
! ==============================================================================

    CASE ( "ncavity" )
      IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization only valid with gas-liq model.')
      END IF ! pRegion%mixtInput%gasModel    
    
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies    
    
      DO icg = 1,pGrid%nCellsTot
        x = pGrid%cofg(XCOORD,icg)

        IF ( x <= 3.0E-03_RFREAL ) THEN
          pCvSpec(1,icg) = 0.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE
          pCvSpec(1,icg) = 1.0_RFREAL 
          pCvSpec(2,icg) = 0.0_RFREAL
        END IF ! x
      END DO ! icg
      
! ==============================================================================
!   Shock-bubble interaction: Quirk and Karni (1996)
! ==============================================================================

    CASE ( "ShockBubble" )
      IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization only valid with gas-liq model.')
      END IF ! pRegion%mixtInput%gasModel   

      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies    

      DO icg = 1,pGrid%nCellsTot
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)
        
        IF ( ((x-0.085_RFREAL)**2 + y**2) <= 6.25E-04_RFREAL ) THEN
          pCvSpec(1,icg) = 0.0_RFREAL
          pCvSpec(2,icg) = 1.0_RFREAL
        ELSE IF ( x < 0.050_RFREAL ) THEN
          pCvSpec(1,icg) = 1.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE
          pCvSpec(1,icg) = 1.0_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        END IF ! x              
      END DO ! icg        
      
! ==============================================================================
!   Shocktube
! ==============================================================================
    
! ------------------------------------------------------------------------------
!   Generic 
! ------------------------------------------------------------------------------
    
    CASE ( "stg1d","stg2d" )
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies
    
      DO icg = 1,pGrid%nCellsTot           
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < pMixtInput%prepRealVal8 ) THEN 
          pCvSpec(1,icg) = pMixtInput%prepRealVal9
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal9
        ELSE IF ( (x >= pMixtInput%prepRealVal8 ) .AND. &
                  (x  < pMixtInput%prepRealVal15) ) THEN 
          pCvSpec(1,icg) = pMixtInput%prepRealVal10
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal10
        ELSE 
          pCvSpec(1,icg) = pMixtInput%prepRealVal16
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal16
        END IF ! x
      END DO ! icg                       
          
! ------------------------------------------------------------------------------
!   Sod 
! ------------------------------------------------------------------------------
    
    CASE ( "st_sod1_mp2","st_sod2_mp2" )
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies
    
      DO icg = 1,pGrid%nCellsTot           
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < 0.0_RFREAL ) THEN 
          pCvSpec(1,icg) = pMixtInput%prepRealVal1
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal1
        ELSE 
          pCvSpec(1,icg) = pMixtInput%prepRealVal2
          pCvSpec(2,icg) = 1.0_RFREAL - pMixtInput%prepRealVal2         
        END IF ! x
      END DO ! icg                       
          
! ==============================================================================
!   Multiphase Riemann problem: Two rarefactions (Toro Case 1) 
! ==============================================================================

    CASE ( "Two_Rarefaction" )
      IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization only valid with gas-liq model.')
      END IF ! pRegion%mixtInput%gasModel   

      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies    

      DO icg = 1,pGrid%nCellsTot
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < 0.5_RFREAL ) THEN
          pCvSpec(1,icg) = 0.01098577_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE
          pCvSpec(1,icg) = 0.01098577_RFREAL
          pCvSpec(2,icg) = 0.0_RFREAL
        END IF ! x
      END DO ! icg        
          
! ==============================================================================
!   Simple volcano model
! ==============================================================================

    CASE ( "volcmod2dv3" )
      IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_TCPERF .OR. & 
           pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_PSEUDO ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization not valid with mixture gas model.')
      END IF ! pRegion%mixtInput%gasModel
    
      DO icg = 1,pGrid%nCellsTot
        y = pGrid%cofg(YCOORD,icg)

        IF ( y < -500.0_RFREAL ) THEN 
          DO iSpec = 1,pRegion%specInput%nSpecies
            pCvSpec(iSpec,icg) = 1.0_RFREAL 
          END DO ! iSpec
        ELSE 
          DO iSpec = 1,pRegion%specInput%nSpecies
            pCvSpec(iSpec,icg) = 0.0_RFREAL
          END DO ! iSpec
        END IF ! y
      END DO ! icg

! ==============================================================================
!   Compressible Shocktube: Air-Air 
! ==============================================================================

    CASE ( "2DShock001" )
      IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Case initialization only valid with gas-liq model.')
      END IF ! pRegion%mixtInput%gasModel   
      
      IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
        WRITE(errorString,'(A,1X,I2)') 'Should be:',pRegion%specInput%nSpecies
        CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__,TRIM(errorString))
      END IF ! pRegion%specInput%nSpecies    

      DO icg = 1,pGrid%nCellsTot
        x = pGrid%cofg(XCOORD,icg)

        IF ( x < 0.5_RFREAL ) THEN
          pCvSpec(1,icg) = 1.0_RFREAL  
          pCvSpec(2,icg) = 0.0_RFREAL
        ELSE
          pCvSpec(1,icg) = 1.0_RFREAL  
          pCvSpec(2,icg) = 0.0_RFREAL
        END IF ! x
      END DO ! icg

! ==============================================================================
!   Default
! ==============================================================================  
  
    CASE DEFAULT     
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%casename

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing species flow field from hard code done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_InitFlowHardCode

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_InitFlowHardCode.F90,v $
! Revision 1.17  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2007/04/05 01:12:43  haselbac
! Added stg1d, modified code to allow 2nd interface
!
! Revision 1.14  2006/05/06 16:51:19  haselbac
! Added kjet2 cases
!
! Revision 1.13  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.12  2006/03/30 20:52:16  haselbac
! Changed ShockBubble hard-code, cosmetics
!
! Revision 1.11  2006/03/26 20:22:18  haselbac
! Added cases for GL model
!
! Revision 1.10  2005/11/17 22:32:32  haselbac
! Added section for gmpjet
!
! Revision 1.9  2005/11/17 14:41:01  haselbac
! Added init for stg2d case
!
! Revision 1.8  2005/11/14 17:02:46  haselbac
! Added support for pseudo-gas model
!
! Revision 1.7  2005/11/10 21:06:03  haselbac
! Changed init for Sod shocktube
!
! Revision 1.6  2005/11/10 02:37:02  haselbac
! Added proper init for Sod shocktube, removed default
!
! Revision 1.5  2005/03/31 17:18:57  haselbac
! Cosmetics only
!
! Revision 1.4  2005/01/30 20:01:05  haselbac
! Added hardcode for Sod shocktube
!
! Revision 1.3  2004/11/12 14:07:37  haselbac
! Modified hard-code for volcano
!
! Revision 1.2  2004/11/10 22:50:29  haselbac
! Added volcano and writing of casename
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************








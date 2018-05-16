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
! Purpose: Set values derived from user input.
!
! Description: None.
!
! Input:
!   regions        Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_DerivedInputValues.F90,v 1.21 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DerivedInputValues(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_DerivedInputValues.F90,v $ $Revision: 1.21 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_DerivedInputValues',&
  'RFLU_DerivedInputValues.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting derived input variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Region-independent values
! ******************************************************************************

  IF ( global%refREnum > 0.0_RFREAL ) THEN
    global%refVisc = global%refDensity*global%RefVelocity*global%refLength/ &
                     global%refREnum
  ELSE
    global%refVisc = 0.0_RFREAL
  END IF ! global%refREnum

! ******************************************************************************
! Region-dependent values
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    pRegion    => regions(iReg)
    pMixtInput => pRegion%mixtInput

! ==============================================================================
!   General
! ==============================================================================

    SELECT CASE ( pMixtInput%fluidModel ) 
      CASE ( FLUID_MODEL_INCOMP ) 
        pMixtInput%nCv = 4
        pMixtInput%nDv = 0 ! TEMPORARY
        pMixtInput%nGv = 0 ! TEMPORARY
      CASE ( FLUID_MODEL_COMP ) 
        pMixtInput%nCv = 5
        pMixtInput%nDv = 3 
        pMixtInput%nGv = 2        

        SELECT CASE ( pMixtInput%gasModel ) 
          CASE ( GAS_MODEL_TCPERF )
            pMixtInput%nGvAct = 0
             
            pMixtInput%indCp  = 0 
            pMixtInput%indMol = 0 
          CASE ( GAS_MODEL_MIXT_TCPERF )
            pMixtInput%nGvAct = 2          
           
            pMixtInput%indCp  = 1 
            pMixtInput%indMol = 1
          CASE ( GAS_MODEL_MIXT_PSEUDO )
            pMixtInput%nGvAct = 2          
           
            pMixtInput%indCp  = 1 
            pMixtInput%indMol = 1 
          CASE ( GAS_MODEL_MIXT_GASLIQ ) 
            pMixtInput%nGvAct = 0

            pMixtInput%indCp  = 1 
            pMixtInput%indMol = 1                     
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pMixtInput%gasModel  
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%fluidModel   
       
    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN 
      pMixtInput%computeTv = .TRUE.
    ELSE 
      pMixtInput%computeTv = .FALSE.
    END IF ! pMixtInput%flowModel

    IF ( pMixtInput%computeTv ) THEN ! transport variables
      pMixtInput%nTv = 2
    ELSE
      pMixtInput%nTv = 0
    END IF ! pMixtInput
        
! ==============================================================================
!   Transport coefficients and dimensionless numbers
! ==============================================================================

    IF ( pMixtInput%refVisc < 0.0_RFREAL ) THEN ! Was not specified
      pMixtInput%refVisc = global%refVisc
    END IF ! pMixtInput%refVisc

    pMixtInput%prLam   = global%prLam
    pMixtInput%prTurb  = global%prTurb
    pMixtInput%scnLam  = global%scnLam
    pMixtInput%scnTurb = global%scnTurb

! ==============================================================================
!   Stencil dimensionality. NOTE if not using one-dimensional stencils, then
!   set stencil dimensionality to be equal to dimensionality.
! ==============================================================================

    IF ( pMixtInput%stencilDimensCells > 1 ) THEN 
      pMixtInput%stencilDimensCells = pMixtInput%dimens
    END IF ! pMixtInput%stencilDimensCells
    
    IF ( pMixtInput%stencilDimensFaces > 1 ) THEN 
      pMixtInput%stencilDimensFaces = pMixtInput%dimens
    END IF ! pMixtInput%stencilDimensFaces
    
    IF ( pMixtInput%stencilDimensBFaces > 1 ) THEN 
      pMixtInput%stencilDimensBFaces = pMixtInput%dimens
    END IF ! pMixtInput%stencilDimensBFaces        

! ==============================================================================
!   Grid speeds
! ==============================================================================

    IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN
      pRegion%grid%indGs = 1
    ELSE
      pRegion%grid%indGs = 0
    END IF ! pMixtInput

! ==============================================================================
!   Time-stepping scheme
! ==============================================================================

    IF ( global%flowType == FLOW_STEADY ) THEN ! steady flow
      global%nrkSteps  = 5
    ELSE ! unsteady flow
      SELECT CASE ( global%rkScheme ) 
        CASE ( RK_SCHEME_4_CLASSICAL ) 
          global%nrkSteps = 4
        CASE ( RK_SCHEME_3_WRAY ) 
          global%nrkSteps = 3
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)            
      END SELECT ! global%rkScheme
    END IF ! global%flowType

    IF ( global%flowType == FLOW_STEADY ) THEN ! steady, upwind scheme
      pMixtInput%ark(1) = 0.0695_RFREAL
      pMixtInput%ark(2) = 0.1602_RFREAL
      pMixtInput%ark(3) = 0.2898_RFREAL
      pMixtInput%ark(4) = 0.5060_RFREAL
      pMixtInput%ark(5) = 1.0000_RFREAL

      pMixtInput%betrk(1) = 1.0000_RFREAL
      pMixtInput%betrk(2) = 0.0000_RFREAL
      pMixtInput%betrk(3) = 0.5600_RFREAL
      pMixtInput%betrk(4) = 0.0000_RFREAL
      pMixtInput%betrk(5) = 0.4400_RFREAL

      pMixtInput%ldiss(1) = 1
      pMixtInput%ldiss(2) = 0
      pMixtInput%ldiss(3) = 1
      pMixtInput%ldiss(4) = 0
      pMixtInput%ldiss(5) = 1
    ELSE ! unsteady
      SELECT CASE ( global%rkScheme ) 
        CASE ( RK_SCHEME_4_CLASSICAL ) 
          pMixtInput%ark(1) = 0.5_RFREAL
          pMixtInput%ark(2) = 0.5_RFREAL
          pMixtInput%ark(3) = 1.0_RFREAL
          pMixtInput%ark(4) = 1.0_RFREAL/6.0_RFREAL

          pMixtInput%betrk(1) = 1.0_RFREAL
          pMixtInput%betrk(2) = 1.0_RFREAL
          pMixtInput%betrk(3) = 1.0_RFREAL      
          pMixtInput%betrk(4) = 1.0_RFREAL            

          pMixtInput%grk(1) = 1.0_RFREAL
          pMixtInput%grk(2) = 2.0_RFREAL
          pMixtInput%grk(3) = 2.0_RFREAL
          pMixtInput%grk(4) = 1.0_RFREAL

          pMixtInput%trk(1) = 0.0_RFREAL
          pMixtInput%trk(2) = 0.5_RFREAL
          pMixtInput%trk(3) = 0.5_RFREAL
          pMixtInput%trk(4) = 1.0_RFREAL

          pMixtInput%ldiss(1) = 1
          pMixtInput%ldiss(2) = 1
          pMixtInput%ldiss(3) = 1
          pMixtInput%ldiss(4) = 1                  

          pMixtInput%cfl = MIN(pMixtInput%cfl,3.0_RFREAL) ! stability margin
        CASE ( RK_SCHEME_3_WRAY ) 
          pMixtInput%ark(1) = 8.0_RFREAL/15.0_RFREAL
          pMixtInput%ark(2) = 5.0_RFREAL/12.0_RFREAL
          pMixtInput%ark(3) = 0.75_RFREAL

          pMixtInput%betrk(1) = 0.0_RFREAL
          pMixtInput%betrk(2) = 0.0_RFREAL
          pMixtInput%betrk(3) = 0.0_RFREAL    

          pMixtInput%grk(1) =  0.0_RFREAL
          pMixtInput%grk(2) = 17.0_RFREAL/25.0_RFREAL
          pMixtInput%grk(3) =  5.0_RFREAL/ 9.0_RFREAL        
 
          pMixtInput%trk(1) = 0.0_RFREAL 
          pMixtInput%trk(2) = 8.0_RFREAL/15.0_RFREAL
          pMixtInput%trk(3) = 2.0_RFREAL/3.0_RFREAL

          pMixtInput%ldiss(1) = 1
          pMixtInput%ldiss(2) = 1
          pMixtInput%ldiss(3) = 1                 

          pMixtInput%cfl = MIN(pMixtInput%cfl,2.0_RFREAL) ! stability margin        
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! global%rkScheme
    END IF ! pMixtInput
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting derived input variables done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DerivedInputValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DerivedInputValues.F90,v $
! Revision 1.21  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.20  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.19  2006/04/13 18:07:09  haselbac
! Bug fix: Changed indCp and indMol for GASLIQ model
!
! Revision 1.18  2006/04/07 14:42:20  haselbac
! Added setting of stencilDimens params
!
! Revision 1.17  2006/03/26 20:21:32  haselbac
! Added support for GL model
!
! Revision 1.16  2006/01/06 22:06:12  haselbac
! Added setting of stencilDimens
!
! Revision 1.15  2005/11/14 16:54:30  haselbac
! Added setting of variables for pseudo-gas model
!
! Revision 1.14  2005/11/10 02:01:30  haselbac
! Added setting of nGvAct, indCp, indMol depending on gasModel
!
! Revision 1.13  2004/11/17 16:27:36  haselbac
! Added setting of RK3 coeffs
!
! Revision 1.12  2004/11/14 19:38:24  haselbac
! Changed nDv and nGv for incompressible fluid model temporarily
!
! Revision 1.11  2004/11/02 02:27:58  haselbac
! Added setting of nCv etc for various fluid models
!
! Revision 1.10  2004/10/19 19:37:30  haselbac
! Cosmetics only
!
! Revision 1.9  2004/04/15 20:18:47  haselbac
! Major bug fix: Accidentally deleted setting of nrkSteps for unsteady flow
!
! Revision 1.8  2004/04/14 02:07:39  haselbac
! Clean-up
!
! Revision 1.7  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/01/29 23:03:32  haselbac
! Removed unnecessary iLev loop and cosmetic changes
!
! Revision 1.5  2003/12/04 03:23:50  haselbac
! Clean-up
!
! Revision 1.4  2003/11/25 21:02:45  haselbac
! Cosmetic changes only
!
! Revision 1.3  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.2  2003/06/19 22:43:59  haselbac
! Changed trk to correspond to different use in rungeKutta
!
! Revision 1.1  2003/01/28 15:53:31  haselbac
! Initial revision, moved from rocflu, use LBOUND and UBOUND
!
! Revision 1.7  2002/10/12 14:56:24  haselbac
! Initialize trk
!
! Revision 1.6  2002/09/09 15:49:00  haselbac
! global now under regions, set nTv for NS, delete nGrad
!
! Revision 1.5  2002/07/25 14:27:43  haselbac
! Added MASTERPROC distinction for output
!
! Revision 1.4  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.3  2002/06/14 20:21:26  haselbac
! Cosmetic changes and prepared addition of grid speed stuff
!
! Revision 1.2  2002/05/04 17:06:31  haselbac
! Minor changes
!
! Revision 1.1  2002/03/26 19:27:52  haselbac
! Initial revision
!
! ******************************************************************************








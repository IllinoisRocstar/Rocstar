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
! Purpose: Compute errors of computed solution relative to exact solution at 
!   probe locations.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************
!
! $Id: RFLU_ComputeExactFlowProbeError.F90,v 1.3 2008/12/06 08:45:05 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeExactFlowProbeError(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
    
  USE RFLU_ModBessel  
  USE RFLU_ModExactFlow, ONLY:  RFLU_ComputeExactFlowPAcoust 
  USE RFLU_ModFlowHardCode, ONLY: RFLU_GetParamsHardCodePAcoust      
    
  USE ModInterfaces, ONLY: MixtPerf_D_CGP, & 
                           MixtPerf_R_CpG
         
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

  LOGICAL :: printErrorNorms
  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  INTEGER :: errorFlag,iBc,icg,im,in,iProbe,iq
  REAL(RFREAL) :: aTot,const,cpGas,dc,de,dInc,dTot,dummyReal,etaqm,gGas, &
                  idc,L,omega,pc,pe,probeTime,pTot,rGas,ro,term,tc,tTot,uc, &
                  ue,vc,ve,wc,we,x,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeExactFlowProbeError.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeExactFlowProbeError', &
                        'RFLU_ComputeExactFlowProbeError.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
      'Computing errors in flow solution at probe locations...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ==============================================================================
! Set pointers
! ==============================================================================

  pGrid      => pRegion%grid
  pCv        => pRegion%mixt%cv
  pDv        => pRegion%mixt%dv
  pMixtInput => pRegion%mixtInput

! ==============================================================================
! Set constants and initialize variables
! ==============================================================================

  cpGas = global%refCp
  gGas  = global%refGamma  
  rGas  = MixtPerf_R_CpG(cpGas,gGas)
  
! ******************************************************************************
! Compute errors in probe quantities 
! ******************************************************************************

  SELECT CASE ( global%casename )

! ==============================================================================
!   Pipe acoustics. NOTE the pipe is assumed to have the x-coordinate 
!   running down the axis. 
! ==============================================================================

    CASE ( "pipeacoust" )
      CALL RFLU_GetParamsHardCodePAcoust(pTot,aTot)
        dTot = MixtPerf_D_CGP(aTot,gGas,pTot)        

        L  = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)) 
        ro = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))        

        im  = MAX(pMixtInput%prepIntVal1,1)
        in  = MAX(pMixtInput%prepIntVal2,1)
        iq  = MAX(pMixtInput%prepIntVal3,1)
        iBc = MAX(MIN(pMixtInput%prepIntVal4,1),0)

        const = MAX(pMixtInput%prepRealVal1,0.0_RFREAL)

        CALL RFLU_JYZOM(im,iq,dummyReal,etaqm,dummyReal,dummyReal)       

        omega = aTot*SQRT((in*global%pi/L)**2 + (etaqm/ro)**2)         

        IF ( global%verbLevel > VERBOSE_LOW ) THEN           
          WRITE(STDOUT,'(A,5X,A,1X,I2)'   ) SOLVER_NAME, &
                'Boundary condition:',iBc          
          WRITE(STDOUT,'(A,5X,A,3(1X,I2))') SOLVER_NAME, &
                'Mode:',im,in,iq                   
          WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                'Total density (kg/m^3):   ',dTot
          WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                'Total pressure (N/m^2):   ',pTot            
          WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                'Angular frequency (rad/s):',omega
          WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                'Constant (-):             ',const                  
        END IF ! global%verbLevel
  
! ------------------------------------------------------------------------------
!       Loop over probes and compute error: NOTE need to rewind because opening 
!       of probe file positions file at end.
! ------------------------------------------------------------------------------
  
        DO iProbe = 1,global%nProbes
          IF ( global%probePos(iProbe,PROBE_REGION) == & 
               pRegion%iRegionGlobal ) THEN 
            icg = global%probePos(iprobe,PROBE_CELL)

            REWIND(IF_PROBE+iProbe-1)

            WRITE(iFileName,'(A,I4.4)') TRIM(global%outDir)// & 
                                        TRIM(global%casename)//'.prbe_',iProbe
            OPEN(IF_EXTR_DATA1,FILE=TRIM(iFileName),FORM='FORMATTED', &
                 STATUS='UNKNOWN',IOSTAT=errorFlag)                            
            global%error = errorFlag
            IF (global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__, &
                             'File: '//TRIM(iFileName))
            END IF ! global%error

            emptyLoop: DO 
              READ(IF_PROBE+iProbe-1,*,IOSTAT=errorFlag) probeTime, & 
                                                         dc,uc,vc,wc,pc,tc
              IF ( errorFlag /= ERR_NONE ) THEN 
                EXIT emptyLoop
              END IF ! errorFlag

              x = pGrid%cofg(XCOORD,icg)
              y = pGrid%cofg(YCOORD,icg)
              z = pGrid%cofg(ZCOORD,icg)

              CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,probeTime,L,ro, & 
                                                iBc,im,in,iq,etaqm,omega, &
                                                dTot,pTot,aTot,const,de,ue, &
                                                ve,we,pe)  
                                                
              WRITE(IF_EXTR_DATA1,'(1PE14.7,3(1X,E13.6))') & 
                    probeTime,pe,pc,((pc-pTot)/(pe-pTot)-1.0_RFREAL)
            END DO emptyLoop
            
            CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)                            
            global%error = errorFlag
            IF (global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__, &
                             'File: '//TRIM(iFileName))
            END IF ! global%error                        
          END IF ! global%probePos  
        END DO ! iProbe

! ==============================================================================
!   Default - due to input error or missing CALL in this routine
! ==============================================================================  

    CASE DEFAULT
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
              'Exact solution not available. Returning to calling procedure.'                                   
      END IF ! global%verbLevel
  END SELECT ! global%casename

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Computing errors in flow solution at probe locations done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeExactFlowProbeError

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeExactFlowProbeError.F90,v $
! Revision 1.3  2008/12/06 08:45:05  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/04/29 12:41:32  haselbac
! Initial revision
!
! ******************************************************************************








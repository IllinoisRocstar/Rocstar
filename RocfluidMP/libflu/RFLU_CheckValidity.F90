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
!   1. Compute and check pressure here because, it being computed later 
!      through a call to MixtureProperties, an invalid value would not be
!      detected.
!
! ******************************************************************************
!
! $Id: RFLU_CheckValidity.F90,v 1.7 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckValidity(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI
  USE ModTools, ONLY: IsNan
  
  USE ModInterfaces, ONLY: MixtPerf_G_CpR, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR, &
                           RFLU_PrintLocInfo
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_INVALID_LOCS = 10
  INTEGER :: icg,indCp,indMol,nLocs
  INTEGER :: loc(MAX_INVALID_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: Eo,gamma,p,rgas,rho,rrho,t,u,v,Vm2,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,gv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckValidity.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckValidity',&
  'RFLU_CheckValidity.F90')

  nLocs = 0

#ifdef ROCPROF
  CALL FPROFILER_BEGINS("RFLU::CheckValidity")
#endif

! ******************************************************************************
! Loop over cells and check for positivity
! ******************************************************************************

  cv => pRegion%mixt%cv
  gv => pRegion%mixt%gv

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol 

  DO icg = 1,pRegion%grid%nCells
    rho  = cv(CV_MIXT_DENS,icg)
    rrho = 1.0_RFREAL/rho
    u    = rrho*cv(CV_MIXT_XMOM,icg)
    v    = rrho*cv(CV_MIXT_YMOM,icg)
    w    = rrho*cv(CV_MIXT_ZMOM,icg)        
    Eo   = rrho*cv(CV_MIXT_ENER,icg)
    rgas = MixtPerf_R_M(gv(GV_MIXT_MOL,icg*indMol))
    gamma= MixtPerf_G_CpR(gv(GV_MIXT_CP,icg*indCp),rgas)
    Vm2  = u*u + v*v + w*w
    p    = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
    t    = MixtPerf_T_DPR(rho,p,rgas)

    IF ( (IsNan(rho) .EQV. .TRUE.) .OR. & 
         (IsNan(u)   .EQV. .TRUE.) .OR. &
         (IsNan(v)   .EQV. .TRUE.) .OR. & 
         (IsNan(w)   .EQV. .TRUE.) .OR. & 
         (IsNan(p)   .EQV. .TRUE.) .OR. &
         (IsNan(t)   .EQV. .TRUE.) ) THEN
      nLocs = nLocs + 1   

      IF ( nLocs == 1 ) THEN 
        WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, & 
              'Invalid variables detected!'
              
        IF ( global%flowType == FLOW_UNSTEADY ) THEN 
          WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                              global%currentTime              
        ELSE 
          WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
                                         'Current iteration number:', &
                                         global%currentIter           
        END IF ! global%flowType                 
                                            
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal 
        WRITE(STDOUT,'(A,6X,A,6(1X,A))') SOLVER_NAME,'#', &
                                         '   Density   ', &
                                         '  x-velocity ', &
                                         '  y-velocity ', &
                                         '  z-velocity ', &
                                         '   Pressure  ', &
                                         ' Temperature '       
      END IF ! nLocs

      IF ( nLocs <= MAX_INVALID_LOCS ) THEN 
        WRITE(STDOUT,'(A,4X,I3,6(1X,E13.6))') SOLVER_NAME,nLocs, & 
                                              rho,u,v,w,p,t
        loc(nLocs,MIN_VAL:MAX_VAL) = icg                                   
      END IF ! nLocs
    END IF ! cv    
  END DO ! icg

! ******************************************************************************
! Write out message and call error handling routine
! ******************************************************************************

  IF ( nLocs > 0 ) THEN 
    IF ( nLocs > MAX_INVALID_LOCS ) THEN 
       WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
             'Only wrote the first',MAX_INVALID_LOCS,'of',nLocs, & 
             'cells with invalid variables.'    
      CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INVALID_LOCS, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    ELSE 
      CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    END IF ! nLocs
    
    CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)   
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF
  CALL FPROFILER_ENDS("RFLU::CheckValidity")
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLU_CheckValidity

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckValidity.F90,v $
! Revision 1.7  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/26 20:21:30  haselbac
! Removed dv declaration and definition
!
! Revision 1.3  2005/07/07 22:44:00  haselbac
! Added profiling calls, cosmetics
!
! Revision 1.2  2004/01/11 02:06:38  jiao
! Eliminated some redundant trailing spaces that made some lines too long.
! This changed was needed to compile with NAG F90 compiler.
!
! Revision 1.1  2003/12/04 03:23:33  haselbac
! Initial revision
!
! ******************************************************************************








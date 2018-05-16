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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfaces.F90,v 1.74 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfaces

  USE ModInterfacesBcond
  USE ModInterfacesIO
  USE ModInterfacesMixt
  USE ModInterfacesUtil

#ifdef RFLO
  USE RFLO_ModInterfacesExternal
  USE RFLO_ModInterfacesLibrary
  USE RFLO_ModInterfacesSolver
#endif

#ifdef RFLU
  USE RFLU_ModInterfacesCommon
#ifdef GENX
  USE RFLU_ModInterfacesExternal
#endif
  USE RFLU_ModInterfacesLibrary
  USE RFLU_ModInterfacesSolver
  USE RFLU_ModInterfacesUtilities
#endif

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE AfterUpdateMP( pRegion,istage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: pRegion
    INTEGER,        INTENT(IN) :: istage
  END SUBROUTINE AfterUpdateMP

  SUBROUTINE AllocateMemoryWork( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE AllocateMemoryWork

  SUBROUTINE BuildVersionString(versionString)
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString
  
  SUBROUTINE CellGradientsMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE CellGradientsMP

  SUBROUTINE CheckPositivityMP( region, pRegion )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE CheckPositivityMP

  SUBROUTINE ConvectiveFluxes( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE ConvectiveFluxes

  SUBROUTINE ConvectiveFluxesMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE ConvectiveFluxesMP

#ifdef RFLO
  SUBROUTINE ExplicitMultistage( regions,ftermNew,residFterm )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    LOGICAL :: ftermNew, residFterm
  END SUBROUTINE ExplicitMultistage
#endif

  SUBROUTINE GlobalCommunicationMP( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE GlobalCommunicationMP

  SUBROUTINE InitCommunicationMP( regions,iReg,istage )
    USE ModDataStruct, ONLY : t_region
    INTEGER, INTENT(IN) :: iReg, istage
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE InitCommunicationMP

  SUBROUTINE IntegrateSourceTermsMP(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE IntegrateSourceTermsMP

  SUBROUTINE NumericalDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE NumericalDissipation

  SUBROUTINE NumericalDissipationMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE NumericalDissipationMP

  SUBROUTINE OptimalLESMP( region, pRegion )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE OptimalLESMP

  SUBROUTINE ReflectPosition(nx,ny,nz,xc,yc,zc,xComp,yComp,zComp)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: nx,ny,nz,xc,yc,zc
    REAL(RFREAL), INTENT(INOUT) :: xComp,yComp,zComp
  END SUBROUTINE ReflectPosition

  SUBROUTINE ReflectVector(nx,ny,nz,xComp,yComp,zComp)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: nx,ny,nz
    REAL(RFREAL), INTENT(INOUT) :: xComp,yComp,zComp
  END SUBROUTINE ReflectVector

  SUBROUTINE RkInitGeneric(region,iStage,icBeg,icEnd,ivBeg,ivEnd,cv,cvOld, &
                           diss)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,diss
    TYPE(t_region) :: region
  END SUBROUTINE RkInitGeneric

  SUBROUTINE RKInitMP( region,istage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
    INTEGER,        INTENT(IN)            :: istage
  END SUBROUTINE RKInitMP

  SUBROUTINE RkInitPointScalar(region,iStage,ivBeg,ivEnd,var,varOld)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: iStage,ivBeg,ivEnd
    REAL(RFREAL), DIMENSION(:), POINTER :: var,varOld
    TYPE(t_region) :: region
  END SUBROUTINE RkInitPointScalar
  
  SUBROUTINE RkInitSD(region,icBeg,icEnd,ivBeg,ivEnd,sd)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icBeg,icEnd,ivBeg,ivEnd
    REAL(RFREAL), DIMENSION(:,:), POINTER :: sd
    TYPE(t_region) :: region
  END SUBROUTINE RkInitSD

  SUBROUTINE RkUpdateGeneric(region,varType,iStage,icBeg,icEnd,ivBeg,ivEnd, & 
                             cv,cvOld,rhs,rhsSum)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
    INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd,varType
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,rhs,rhsSum
  END SUBROUTINE RkUpdateGeneric

  SUBROUTINE RKUpdateMP( region, iReg, istage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
    INTEGER,        INTENT(IN)            :: iReg, istage
  END SUBROUTINE RKUpdateMP

  SUBROUTINE RkUpdatePointScalar(region,iStage,ivBeg,ivEnd,var,varOld, &
                                 rhs,rhsSum)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
    INTEGER, INTENT(IN) :: iStage,ivBeg,ivEnd
    REAL(RFREAL), DIMENSION(:), POINTER :: var,varOld,rhs,rhsSum
  END SUBROUTINE RkUpdatePointScalar

  SUBROUTINE RungeKuttaMP( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RungeKuttaMP

  SUBROUTINE SourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE SourceTerms

  SUBROUTINE SourceTermsMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE SourceTermsMP

  SUBROUTINE UpdateDependentVarsMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE UpdateDependentVarsMP

  SUBROUTINE ViscousFluxes( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE ViscousFluxes

  SUBROUTINE ViscousFluxesMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE ViscousFluxesMP

  SUBROUTINE ZeroResidualsMP(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE ZeroResidualsMP

  END INTERFACE

END MODULE ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfaces.F90,v $
! Revision 1.74  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.73  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.72  2006/08/19 15:38:53  mparmar
! Added RkInitPointScalar and RkUpdatePointScalar
!
! Revision 1.71  2006/01/06 22:08:56  haselbac
! Enclosed if for ExplicitMultistage within ifdef RFLO
!
! Revision 1.70  2004/12/02 15:26:05  haselbac
! Added entry for BuildVersionString
!
! Revision 1.69  2004/11/29 17:38:28  wasistho
! get rid of physical modules ModInterfaces
!
! Revision 1.68  2004/11/29 17:13:02  wasistho
! commented physical module ModInterfaces
!
! Revision 1.67  2004/11/17 16:29:23  haselbac
! Adapted interface for rkUpdateGeneric
!
! Revision 1.66  2004/08/02 23:12:52  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.65  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.64  2004/05/05 20:49:43  fnajjar
! Added entry for ReflectVector and corrected ReflectPosition
!
! Revision 1.63  2004/04/08 01:34:54  haselbac
! Added interface for ReflectVector
!
! Revision 1.62  2004/04/01 21:28:15  haselbac
! Added entry for integrateSourceTermsMP
!
! Revision 1.61  2004/03/25 21:14:21  jferry
! changed AfterUpdate to call most subroutines only after final RK stage
!
! Revision 1.60  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.59  2004/02/26 21:14:58  wasistho
! added globalCommunication
!
! Revision 1.58  2004/01/31 03:57:57  haselbac
! Removed interface for rungeKutta, no longer needed
!
! Revision 1.57  2004/01/29 22:57:23  haselbac
! Added interface for updateDependentVarsMP.F90
!
! Revision 1.56  2003/12/04 03:28:26  haselbac
! Added interface for CellGradientsMP
!
! Revision 1.55  2003/11/25 21:03:11  haselbac
! Added interfaces for generic MP routines
!
! Revision 1.54  2003/09/26 21:43:23  fnajjar
! Commented ModInterfaces pertinent to Physical modules
!
! Revision 1.53  2003/08/28 20:29:38  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.52  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
! Revision 1.51  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
!******************************************************************************







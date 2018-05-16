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
! Purpose: Obtain viscous fluxes based the Bardina scale similarity model. 
!
! Description: Get the scale-similarity tauij then add the viscous fluxes to 
!              the dissipation residual, mixt%diss.
!              The ss tauij is obtained by calling LesLij then
!              the viscous fluxes by calling VFluxHybrid.
!
! Input: region  = data of current region
!        ibn,ien = begin and end node index
!
! Output: tauij, eddy viscosity (mueT, derived from tauij) and viscous fluxes.
!
! Notes: tauij in i, j and k faces are kept in fISij, fJSij and fKSij.
!
!******************************************************************************
!
! $Id: TURB_LesFluxScalSim.F90,v 1.9 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesFluxScalSim( region,ibn,ien )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModTurbulence, ONLY : t_turb
  USE ModGlobal, ONLY     : t_global
#ifdef RFLU
  USE ModBndPatch, ONLY   : t_patch
  USE TURB_ModInterfaces, ONLY : TURB_FluLesBLij, TURB_FluLesC2F
#endif
#ifdef RFLO
  USE TURB_ModInterfaces, ONLY : TURB_FloLesGenC2F
#endif
  USE TURB_ModInterfaces, ONLY : TURB_VFluxHybrid, TURB_LesLij
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif
  INTEGER                 :: ibn, ien

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
#ifdef RFLU
  TYPE(t_patch), POINTER  :: patch
#endif

  INTEGER :: errorFlag
#ifdef RFLO              
  INTEGER :: iLev
#endif
#ifdef RFLU              
  INTEGER :: nPatches, nBFaces
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesFluxScalSim.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesFluxScalSim',&
  'TURB_LesFluxScalSim.F90' )

! get pointers and parameters ------------------------------------------------
#ifdef RFLO
  iLev =  region%currLevel
  turb => region%levels(iLev)%turb
#endif
#ifdef RFLU
  nPatches = region%grid%nPatches
  nBFaces  = 0

  DO iPatch = 1,nPatches
    patch   => region%patches(iPatch)
    nBFaces =  nBFaces + patch%nBTris + patch%nBQuads
  END DO ! iPatch
  turb => region%turb
#endif

! get viscous fluxes from scale similarity contribution; first allocate arrays

  ALLOCATE( turb%fVar(CV_TURB_NELM,ibn:ien),stat=errorFlag )    
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%ffVar(CV_TURB_NELM,ibn:ien),stat=errorFlag )    
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%fISij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

#ifdef RFLO
  ALLOCATE( turb%fJSij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%fKSij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif
#ifdef RFLU
  ALLOCATE( turb%bfVar(CV_TURB_NELM,nBFaces),stat=errorFlag )    
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%bffVar(CV_TURB_NELM,nBFaces),stat=errorFlag )    
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%bfISij(TENSOR_SYMM_NELM,nBFaces),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif

! Leonard stress lij in I direction
#ifdef RFLO
  CALL TURB_FloLesGenC2F( region,DIRI )
#endif
#ifdef RFLU
  CALL TURB_FluLesC2F( region )
#endif

  CALL TURB_LesLij( region,DIRI,region%turbInput%filterWidth,turb%fISij )
#ifdef RFLU
  CALL TURB_FluLesBLij( region,region%turbInput%filterWidth,turb%bfISij )
#endif

#ifdef RFLO
! Leonard stress lij in J direction
  CALL TURB_FloLesGenC2F( region,DIRJ )
  CALL TURB_LesLij( region,DIRJ,region%turbInput%filterWidth,turb%fJSij )

! Leonard stress lij in K direction
  CALL TURB_FloLesGenC2F( region,DIRK )
  CALL TURB_LesLij( region,DIRK,region%turbInput%filterWidth,turb%fKSij )
#endif

#ifdef RFLU
!  apply boundary conditions for LES variables
!  CALL TURB_FluLesBndConditions( targetFlag )
#endif

! get viscous fluxes
  CALL TURB_VFluxHybrid( region )

! deallocate retired arrays
  DEALLOCATE( turb%fVar,turb%ffVar,turb%fISij )
#ifdef RFLO
  DEALLOCATE( turb%fJSij,turb%fKSij )
#endif
#ifdef RFLU
  DEALLOCATE( turb%bfVar,turb%bffVar,turb%bfISij )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesFluxScalSim

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesFluxScalSim.F90,v $
! Revision 1.9  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2004/07/30 22:34:45  wasistho
! replaced floLesUniC2F by floLesGenC2F
!
! Revision 1.6  2004/05/28 01:59:12  wasistho
! update unstructured grid LES
!
! Revision 1.5  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.4  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/19 02:50:09  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.3  2003/08/29 01:41:06  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.2  2002/10/16 01:59:49  wasistho
! Changed global%error flag
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************








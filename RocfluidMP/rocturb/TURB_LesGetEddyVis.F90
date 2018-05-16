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
! Purpose: Get eddy viscosity based the dynamic LES model.
!
! Description: First, cell filtered cv is computed, then the filtered 
!              velocities are obtained from this. The gradients of the 
!              filtered velocities are computed. This is in turn used
!              to construct the strain rate tensor of the filtered field.
!              The strain rate tensor of filtered field is employed in the 
!              formulation of the eddy viscosity of the dynamic models.
!
! Input: region  = data of current region
!        ibn,ien = begin and end node index 
!
! Output: eddy viscosity (mueT), stored in turb%mueT(DIRI:DIRK,:)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesGetEddyVis.F90,v 1.13 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesGetEddyVis( region,ibc,iec,ibn,ien )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModTurbulence, ONLY : t_turb
  USE ModGlobal, ONLY     : t_global
#ifdef RFLU
  USE ModBndPatch, ONLY   : t_patch
  USE RFLU_ModDifferentiationFaces, ONLY: RFLU_ComputeGradFacesWrapper
  USE RFLU_ModDifferentiationBFaces, ONLY: RFLU_ComputeGradBFacesWrapper
  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace
#endif
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_CalcGradVector
#endif
  USE TURB_ModInterfaces, ONLY : TURB_LesTestRhoV, TURB_CalcStrainRate, &
                                 TURB_LesCalcEddyVis
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region),  TARGET :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif
  INTEGER :: ibc, iec, ibn, ien

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
#ifdef RFLU
  TYPE(t_patch), POINTER  :: pPatch
#endif
  
  INTEGER               :: iBegV, iEndV, iBegG, iEndG, nGrad, errorFlag
  REAL(RFREAL), POINTER :: ccVar(:,:)
#ifdef RFLO
  INTEGER               :: iLev, gradIndx(9)
  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
#endif
#ifdef RFLU
  INTEGER               :: nPatches, nBFaces, nBFacesTot, gradIndx(3)
  REAL(RFREAL), POINTER :: gradi(:,:,:), bGradi(:,:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesGetEddyVis.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesGetEddyVis',&
  'TURB_LesGetEddyVis.F90' )

! get pointers and parameters ------------------------------------------------
  nGrad =  region%turbInput%nGrad
#ifdef RFLO
  iLev  =  region%currLevel
  turb  => region%levels(iLev)%turb
#endif
#ifdef RFLU
  nPatches   = region%grid%nPatches
  nBFaces    = 0
  nBFacesTot = 0

  DO iPatch = 1,nPatches
    pPatch => region%patches(iPatch)

    nBFaces    = nBFaces    + pPatch%nBTris    + pPatch%nBQuads
    nBFacesTot = nBFacesTot + pPatch%nBTrisTot + pPatch%nBQuadsTot
  END DO ! iPatch
  turb  => region%turb
#endif

! allocate and point arrays required by LesTestRhoV and CalcGradVector

  ALLOCATE( turb%ccVar(CV_TURB_NELM,ibc:iec),stat=errorFlag ) 
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

#ifdef RFLO
  ALLOCATE( turb%gradi(nGrad,ibn:ien),stat=errorFlag )         
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%gradj(nGrad,ibn:ien),stat=errorFlag )         
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%gradk(nGrad,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif
#ifdef RFLU
  ALLOCATE( turb%gradi(XCOORD:ZCOORD,nGrad,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%bGradi(XCOORD:ZCOORD,nGrad,nBFaces),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif

  ccVar  => turb%ccVar
  gradi  => turb%gradi
#ifdef RFLO
  gradj  => turb%gradj
  gradk  => turb%gradk
#endif
#ifdef RFLU
  bGradi => turb%bGradi
#endif

! compute cell-filtered cv and store in ccVar
  CALL TURB_LesTestRhoV( region,ibc,iec )

! compute gradients of ccVar
#ifdef RFLO
  iBegV =  CV_TURB_UVEL
  iEndV =  CV_TURB_WVEL
  iBegG =  GR_TURB_UX
  iEndG =  GR_TURB_WZ
  CALL RFLO_CalcGradVector( region,iBegV,iEndV,iBegG,iEndG,ccVar, &
                            gradi,gradj,gradk )
#endif
#ifdef RFLU
  iBegV  =  CV_TURB_UVEL
  iEndV  =  CV_TURB_WVEL
  iBegG  =  GR_TURB_UX
  iEndG  =  GR_TURB_WX
  CALL RFLU_ComputeGradFacesWrapper( region,iBegV,iEndV,iBegG,iEndG,ccVar,gradi )
  DO iPatch = 1,nPatches
    pPatch => region%patches(iPatch)

    IF ( RFLU_DecideNeedBGradFace(region,pPatch) .EQV. .TRUE. ) THEN
      CALL RFLU_ComputeGradBFacesWrapper( region,pPatch,iBegV,iEndV,iBegG,iEndG, &
                                          ccVar,bGradi )
    END IF ! RFLU_DecideNeedBGradFace
  ENDDO ! iPatch 
#endif
  DEALLOCATE( turb%ccVar )

! compute sij based on cell-filtered cv, Sij(v) = Sij(test[bar(u_i)]) 
! and store this in fISij, fJSij, fKSij
  gradIndx(1) = GR_TURB_UX
  gradIndx(2) = GR_TURB_VX
  gradIndx(3) = GR_TURB_WX
#ifdef RFLO
  gradIndx(4) = GR_TURB_UY
  gradIndx(5) = GR_TURB_VY
  gradIndx(6) = GR_TURB_WY
  gradIndx(7) = GR_TURB_UZ
  gradIndx(8) = GR_TURB_VZ
  gradIndx(9) = GR_TURB_WZ

  CALL TURB_CalcStrainRate( region,ibn,ien,gradIndx,gradi,gradj,gradk, &
                            turb%fISij,turb%fJSij,turb%fKSij )
  DEALLOCATE( turb%gradi,turb%gradj,turb%gradk )
#endif
#ifdef RFLU
  CALL TURB_CalcStrainRate( region, ibn,ien, gradIndx, gradi, turb%fISij )
  IF (nPatches > 0) &
  CALL TURB_CalcStrainRate( region,1,nBFaces,gradIndx,bGradi,turb%bfISij )
  DEALLOCATE( turb%gradi, turb%bGradi )
#endif

! obtain eddyviscosity in i, j and k faces
  CALL TURB_LesCalcEddyVis( region,ibn,ien,DIRI )
#ifdef RFLO
  CALL TURB_LesCalcEddyVis( region,ibn,ien,DIRJ )
  CALL TURB_LesCalcEddyVis( region,ibn,ien,DIRK )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesGetEddyVis

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesGetEddyVis.F90,v $
! Revision 1.13  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2006/08/19 15:40:36  mparmar
! Added use of RFLU_DecideNeedBGradFace
!
! Revision 1.10  2006/04/07 15:06:06  haselbac
! Bug fix: Incorrect ifs
!
! Revision 1.9  2006/04/07 14:56:02  haselbac
! Adapted to changes in f and bf grad routines
!
! Revision 1.8  2005/12/20 20:43:29  wasistho
! adapted to changing in Rocflu on face gradient routines
!
! Revision 1.7  2004/05/28 01:58:49  wasistho
! update unstructured grid LES
!
! Revision 1.6  2004/05/17 20:33:39  wasistho
! reordering gradIndx for more efficient cache memory addressing
!
! Revision 1.5  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.4  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/19 02:50:18  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.5  2003/10/09 23:07:41  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.4  2003/08/29 01:42:25  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.3  2003/05/16 05:43:34  wasistho
! modified array range of CC-filtered
!
! Revision 1.2  2002/10/16 02:00:06  wasistho
! Changed global%error flag
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************








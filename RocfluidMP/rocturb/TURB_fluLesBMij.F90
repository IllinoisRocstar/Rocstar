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
! Purpose: Compute the mij stress appearing in the dynamic LES models.
!
! Description: We determine mij as appearing in the dynamic formulation:
!              M_{ij}=-test[bar(rho)](kappa \Delta)^2 |S(v)|S_{ij(v)
!                     +test[bar(rho)\Delta^{2}|S(u)|S_{ij}(u)]
!              in which v=test[bar(rho ui)]/test[bar(rho)] where test denotes
!              the test-filter with twice the filter-width. M_ij is computed
!              on the cell faces. S_ij(u) and Sij(test[u]) are already 
!              available, respectively, in mSij for mixture field and fSij for 
!              filtered field at cell faces. Moreover, we use the fact that we 
!              know 1/test[bar(rho)] (rhoFTb) and 1/[bar(rho)] (rhoF) at faces.
!
! Input: region = data of current region
!        ijk    = ijk-face is being treated
!
! Output: mij = resulting stress tensor in dynamic LES models
!
! Notes: This routine relevant only for physical boundaries of FLU regions.
!
!******************************************************************************
!
! $Id: TURB_fluLesBMij.F90,v 1.3 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FluLesBMij( region,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
!  USE RFLU_ModInterpolation, ONLY : RFLU_InterpFaces2FacesPatches
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: region
  INTEGER :: ijk

! ... loop variables
  INTEGER :: i, j, k, ijkN, ijkC0, ijkC1, iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

  INTEGER      :: ibn,ien,idBeg,idEnd,nPatches,nBFaces,tNDel(DIRI:DIRK)
  REAL(RFREAL) :: kappa2,rkappa2,two3rd
  REAL(RFREAL) :: delFac2,delFKap2,forTerm
  REAL(RFREAL), POINTER :: bfVol(:),bfSij(:,:),bfVar(:,:),bffVar(:,:)
  REAL(RFREAL), POINTER :: bSij(:,:),bMij(:,:)
  REAL(RFREAL), POINTER :: field(:,:),tField(:,:)
  REAL(RFREAL), ALLOCATABLE :: modStrain(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_fluLesBMij.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FluLesBMij',&
  'TURB_fluLesBMij.F90' )

! get indices and pointers ----------------------------------------------------

  nPatches = region%grid%nPatches
  nBFaces  = 0

  DO iPatch = 1,nPatches
    patch   => region%patches(iPatch)
    nBFaces =  nBFaces + patch%nBTris + patch%nBQuads
  END DO ! iPatch
  ibn = 1
  ien = nBFaces

  bfVar  => region%turb%bfVar
  bffVar => region%turb%bffVar
  bMij   => region%turb%bMij
  bfVol  => region%turb%bfVolI
  bSij   => region%turb%bmIsij
  bfSij  => region%turb%bfISij

! allocate work arrays used within this scope

  ALLOCATE( modStrain(ibn:ien) )
  ALLOCATE( field(E11:E33,ibn:ien),tField(E11:E33,ibn:ien) )

! test filter width is twice bar filter width

  tNDel(DIRI) = 2*region%turbInput%filterWidth(DIRI)

! filterwidth on the fg-level (test-bar) is sqrt(5) times larger than 
! on the f-level (bar), i.e. kappa=sqrt(5)

  kappa2 = 5._RFREAL
  rkappa2= 1._RFREAL/kappa2
  two3rd = 2._RFREAL/3._RFREAL

! determine the bar filter width delta which appears in the dynamic model; 
! this factor in square is defined as (kappa*delfac)^2 and stored in delFKap2 

  delFac2  = region%turbInput%delFac2
  delFKap2 = delFac2*kappa2

! first calculate test[bar(rho)](\kappa \Delta)^2|S(v)|Sij(v) at cell faces; 
! we obtain Sij(v) from array fSij; 

  DO ijkN=ibn,ien

! - generate the 6-th component of fSij using the tracelessness of S_{ij} 

    bfSij(E33,ijkN)=-(bfSij(E11,ijkN)+bfSij(E22,ijkN))

! - get modulus of Sij(v): |S(v)|
    modStrain(ijkN)=sqrt(bfSij(E11,ijkN)*bfSij(E11,ijkN) &
                        +bfSij(E22,ijkN)*bfSij(E22,ijkN) &
                        +bfSij(E11,ijkN)*bfSij(E22,ijkN) &
                        +bfSij(E12,ijkN)*bfSij(E12,ijkN) &
                        +bfSij(E13,ijkN)*bfSij(E13,ijkN) &
                        +bfSij(E23,ijkN)*bfSij(E23,ijkN))

! - in turn, (kappa*delta)^2 in fg-level is stored in tField

    tField(E11,ijkN) = delFKap2*bfVol(ijkN)**two3rd

! - fg-level contribution to M_{ij} can now be computed since
!   we already have test(bar[rho]) available at cell faces,
!   in fact we have 1/test(bar[rho]), so we devide instead of multiply;
!   install the first term of mij; the space of fSij in used in the process

    forTerm       = -tField(E11,ijkN)/bffVar(CV_TURB_DENS,ijkN)* &
                     modStrain(ijkN)

    bMij(E11,ijkN) = forTerm*bfSij(E11,ijkN)
    bMij(E12,ijkN) = forTerm*bfSij(E12,ijkN)
    bMij(E13,ijkN) = forTerm*bfSij(E13,ijkN)
    bMij(E22,ijkN) = forTerm*bfSij(E22,ijkN)
    bMij(E23,ijkN) = forTerm*bfSij(E23,ijkN)
    bMij(E33,ijkN) = forTerm*bfSij(E33,ijkN)

! - as final step we determine contribution of test-filtered f-level term;
!   note we have bar[rho] available at cell faces bRho stored in fVar; 
!   in fact we have 1/bRho in fVar(CV_TURB_DENS,:), so we devide by bRho
!   instead of multiply and compute modulus of Sij: |S(u)|

    modStrain(ijkN) = sqrt(bSij(E11,ijkN)*bSij(E11,ijkN) &
                          +bSij(E22,ijkN)*bSij(E22,ijkN) &
                          +bSij(E11,ijkN)*bSij(E22,ijkN) &
                          +bSij(E12,ijkN)*bSij(E12,ijkN) &
                          +bSij(E13,ijkN)*bSij(E13,ijkN) &
                          +bSij(E23,ijkN)*bSij(E23,ijkN))

! - finish computation of mij; filterwidth^2 at f-level is tField/kappa2;
!   build terms need to be filtered

    forTerm = tField(E11,ijkN)/bfVar(CV_TURB_DENS,ijkN)*rkappa2* &
              modStrain(ijkN)

    field(E11,ijkN) =  forTerm*bSij(E11,ijkN)
    field(E12,ijkN) =  forTerm*bSij(E12,ijkN)
    field(E13,ijkN) =  forTerm*bSij(E13,ijkN)
    field(E22,ijkN) =  forTerm*bSij(E22,ijkN)
    field(E23,ijkN) =  forTerm*bSij(E23,ijkN)
    field(E33,ijkN) = -forTerm*(bSij(E11,ijkN)+bSij(E22,ijkN))
  ENDDO

! perform test filtering:

  idBeg = E11
  idEnd = E33

!  CALL RFLU_InterpFaces2FacesPatches( region,field,tField )

! and complete components of mij:

  DO ijkN=ibn,ien
    bMij(E11,ijkN) = bMij(E11,ijkN) + tField(E11,ijkN)
    bMij(E12,ijkN) = bMij(E12,ijkN) + tField(E12,ijkN)
    bMij(E13,ijkN) = bMij(E13,ijkN) + tField(E13,ijkN)
    bMij(E22,ijkN) = bMij(E22,ijkN) + tField(E22,ijkN)
    bMij(E23,ijkN) = bMij(E23,ijkN) + tField(E23,ijkN)
    bMij(E33,ijkN) = bMij(E33,ijkN) + tField(E33,ijkN)
  ENDDO

! deallocate temporary arrays

  DEALLOCATE( modStrain,field,tField )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FluLesBMij

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_fluLesBMij.F90,v $
! Revision 1.3  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/05/28 02:03:58  wasistho
! update unstructured grid LES
!
!
!
!******************************************************************************








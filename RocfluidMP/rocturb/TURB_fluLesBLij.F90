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
! Purpose: Compute the scale similarity turbulent stress tensor, which
!          is the same as the Leonard stress term in the dynamic LES
!          models. SS-tauij uses bar-filter while Leonard stress uses
!          test-filter.
!
! Description: The turbulent stress tensor follows the SS model and can
!              be seen as Leonard tensor if the outer bar filter is replaced
!              by test filter. This tensor is defined by
!              l_{ij}=bar(bar[rho*ui]bar[rho*uj]/bar[rho])
!                    -bar(bar[rho*ui])bar(bar[rho*uj])/bar(bar[rho])
!              where bar denotes the bar-filter with filter-width defined by
!              nDel(DIRI:DIRK). On entry we have the solution components 
!              at cell faces and perform the filtering from face to face. 
!              We return the resulting tensor at all cell faces including  
!              dummies. Also 1/bar(bar[rho]) is available stored ffVar(1,:).
!
! Input: region = data of current region
!        nDel   = three components filter width paramater
!
! Output: lij = resulting stress tensor of scale similarity or Leonard term
!
! Notes: This routine relevant only for physical boundaries of FLU regions.
!  1. Removed CALL to RFLU_InterpCells2FacesPatches because it no longer exists 
!     This routine might become broken because of this
!
!******************************************************************************
!
! $Id: TURB_fluLesBLij.F90,v 1.5 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FluLesBLij( region,nDel,lij )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: region
  INTEGER                 :: nDel(DIRI:DIRK)
  REAL(RFREAL), POINTER   :: lij(:,:)

! ... loop variables
  INTEGER :: ijkN, iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

  INTEGER :: ibn,ien,idBeg,idEnd,nPatches,nBFaces
  REAL(RFREAL), POINTER :: bfVar(:,:),bffVar(:,:)
  REAL(RFREAL), POINTER :: tens(:,:),tensBar(:,:)
  REAL(RFREAL), POINTER :: cv(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_fluLesBLij.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FluLesBLij',&
  'TURB_fluLesBLij.F90' )

! get indices and pointers --------------------------------------------------

  nPatches = region%grid%nPatches
  nBFaces  = 0

  DO iPatch = 1,nPatches
    patch   => region%patches(iPatch)
    nBFaces =  nBFaces + patch%nBTris + patch%nBQuads
  END DO ! iPatch
  ibn = 1
  ien = nBFaces

  cv     => region%mixt%cv
  bfVar  => region%turb%bfVar
  bffVar => region%turb%bffVar

! allocate temporary arrays

  ALLOCATE( tens(E11:E33,ibn:ien),tensBar(E11:E33,ibn:ien) )

! ff-filtering bRhoF and bRhoUiF (i=1,3) stored in fVar to get 
! bar[bRhoF] and bar[bRhoUiF] stored in fbVar

  idBeg = CV_TURB_DENS
  idEnd = CV_TURB_ZMOM

! for efficiency we store 1/bRhoF in bRhoF and 1/bbRhoF in bbRhoF;
! next, build components of scale-similarity (Bardina) turbulent stress 
! tensor lij-component: store bRhoUi*bRhoUj/bRho in tens

  DO ijkN = ibn,ien
    bfVar(CV_TURB_DENS,ijkN)=1._RFREAL/bfVar(CV_TURB_DENS,ijkN)
    bffVar(CV_TURB_DENS,ijkN)=1._RFREAL/bffVar(CV_TURB_DENS,ijkN)

    tens(E11,ijkN) = bfVar(CV_TURB_XMOM,ijkN)*bfVar(CV_TURB_XMOM,ijkN)* &
                     bfVar(CV_TURB_DENS,ijkN)
    tens(E12,ijkN) = bfVar(CV_TURB_XMOM,ijkN)*bfVar(CV_TURB_YMOM,ijkN)* &
                     bfVar(CV_TURB_DENS,ijkN)
    tens(E13,ijkN) = bfVar(CV_TURB_XMOM,ijkN)*bfVar(CV_TURB_ZMOM,ijkN)* &
                     bfVar(CV_TURB_DENS,ijkN)
    tens(E22,ijkN) = bfVar(CV_TURB_YMOM,ijkN)*bfVar(CV_TURB_YMOM,ijkN)* &
                     bfVar(CV_TURB_DENS,ijkN)
    tens(E23,ijkN) = bfVar(CV_TURB_YMOM,ijkN)*bfVar(CV_TURB_ZMOM,ijkN)* &
                     bfVar(CV_TURB_DENS,ijkN)
    tens(E33,ijkN) = bfVar(CV_TURB_ZMOM,ijkN)*bfVar(CV_TURB_ZMOM,ijkN)* &
                     bfVar(CV_TURB_DENS,ijkN)
  ENDDO

! ff-filter bRhoUi*bRhoUj/bRho stored in tens to get 
! bar[bRhoUi*bRhoUj/bRho] in tensBar

  idBeg = E11
  idEnd = E33

!  CALL RFLU_InterpFaces2FacesPatches( region,tens,tensBar )

! combine and find lij=bar[bRhoUi*bRhoUj/bRho]-bbRhoUi*bbRhoUj/bbRho

  DO ijkN=ibn,ien
    lij(E11,ijkN)=tensBar(E11,ijkN)-bffVar(CV_TURB_DENS,ijkN)* &
                  bffVar(CV_TURB_XMOM,ijkN)*bffVar(CV_TURB_XMOM,ijkN)

    lij(E12,ijkN)=tensBar(E12,ijkN)-bffVar(CV_TURB_DENS,ijkN)* &
                  bffVar(CV_TURB_XMOM,ijkN)*bffVar(CV_TURB_YMOM,ijkN)

    lij(E13,ijkN)=tensBar(E13,ijkN)-bffVar(CV_TURB_DENS,ijkN)* &
                  bffVar(CV_TURB_XMOM,ijkN)*bffVar(CV_TURB_ZMOM,ijkN)

    lij(E22,ijkN)=tensBar(E22,ijkN)-bffVar(CV_TURB_DENS,ijkN)* &
                  bffVar(CV_TURB_YMOM,ijkN)*bffVar(CV_TURB_YMOM,ijkN)

    lij(E23,ijkN)=tensBar(E23,ijkN)-bffVar(CV_TURB_DENS,ijkN)* &
                  bffVar(CV_TURB_YMOM,ijkN)*bffVar(CV_TURB_ZMOM,ijkN)

    lij(E33,ijkN)=tensBar(E33,ijkN)-bffVar(CV_TURB_DENS,ijkN)* &
                  bffVar(CV_TURB_ZMOM,ijkN)*bffVar(CV_TURB_ZMOM,ijkN)
  ENDDO

! deallocate temporary arrays

  DEALLOCATE( tens,tensBar )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FluLesBLij

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_fluLesBLij.F90,v $
! Revision 1.5  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:41:03  mparmar
! Removed call to RFLU_InterpCells2FacesPatches
!
! Revision 1.2  2004/08/07 01:07:18  wasistho
! bugfixed, defined pointer cv
!
! Revision 1.1  2004/05/28 02:03:58  wasistho
! update unstructured grid LES
!
!
!
!******************************************************************************








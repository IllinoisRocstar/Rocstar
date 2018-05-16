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
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesMij.F90,v 1.8 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesMij( region,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLU
!  USE RFLU_ModInterpolation, ONLY: RFLU_InterpFaces2Faces
#endif
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesUniFiltFF, TURB_FloLesGenFiltFF

#include "Indexing.h"
#endif
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region)          :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif
  INTEGER :: ijk

! ... loop variables
  INTEGER :: i, j, k, ijkN, ijkC0, ijkC1

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER      :: ibn,ien,idBeg,idEnd,tNDel(DIRI:DIRK)
  REAL(RFREAL) :: kappa2,rkappa2,two3rd
  REAL(RFREAL) :: delFac2,delFKap2,forTerm
  REAL(RFREAL), POINTER :: fvol(:),fSij(:,:),fVar(:,:),ffVar(:,:)
  REAL(RFREAL), POINTER :: sij(:,:),mij(:,:)
  REAL(RFREAL), POINTER :: field(:,:),tField(:,:)
  REAL(RFREAL), ALLOCATABLE :: modStrain(:)

#ifdef RFLO
  INTEGER      :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER      :: iLev,iNOff,ijNOff
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesMij.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesMij',&
  'TURB_LesMij.F90' )

! get indices and pointers --------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  fVar  => region%levels(iLev)%turb%fVar
  ffVar => region%levels(iLev)%turb%ffVar
  mij   => region%levels(iLev)%turb%mij
  ibn   =  IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien   =  IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  IF (ijk==DIRI) THEN
    fvol => region%levels(iLev)%turb%fvolI
    sij  => region%levels(iLev)%turb%mIsij
    fSij => region%levels(iLev)%turb%fISij
  ELSEIF (ijk==DIRJ) THEN
    fvol => region%levels(iLev)%turb%fvolJ
    sij  => region%levels(iLev)%turb%mJsij
    fSij => region%levels(iLev)%turb%fJSij
  ELSEIF (ijk==DIRK) THEN
    fvol => region%levels(iLev)%turb%fvolK
    sij  => region%levels(iLev)%turb%mKsij
    fSij => region%levels(iLev)%turb%fKSij
  ENDIF
#endif
#ifdef RFLU
  fVar  => region%turb%fVar
  ffVar => region%turb%ffVar
  mij   => region%turb%mij
  fvol  => region%turb%fvolI
  sij   => region%turb%mIsij
  fSij  => region%turb%fISij
  ibn   =  1
  ien   =  region%grid%nFaces
#endif

! allocate work arrays used within this scope

  ALLOCATE( modStrain(ibn:ien) )
  ALLOCATE( field(E11:E33,ibn:ien),tField(E11:E33,ibn:ien) )

! test filter width is twice bar filter width

  tNDel(DIRI) = 2*region%turbInput%filterWidth(DIRI)
#ifdef RFLO
  tNDel(DIRJ) = 2*region%turbInput%filterWidth(DIRJ)
  tNDel(DIRK) = 2*region%turbInput%filterWidth(DIRK)
#endif

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

    fSij(E33,ijkN)=-(fSij(E11,ijkN)+fSij(E22,ijkN))

! - get modulus of Sij(v): |S(v)|
    modStrain(ijkN)=sqrt(fSij(E11,ijkN)*fSij(E11,ijkN) &
                        +fSij(E22,ijkN)*fSij(E22,ijkN) &
                        +fSij(E11,ijkN)*fSij(E22,ijkN) &
                        +fSij(E12,ijkN)*fSij(E12,ijkN) &
                        +fSij(E13,ijkN)*fSij(E13,ijkN) &
                        +fSij(E23,ijkN)*fSij(E23,ijkN))

! - in turn, (kappa*delta)^2 in fg-level is stored in tField

    tField(E11,ijkN) = delFKap2*fvol(ijkN)**two3rd

! - fg-level contribution to M_{ij} can now be computed since
!   we already have test(bar[rho]) available at cell faces,
!   in fact we have 1/test(bar[rho]), so we devide instead of multiply;
!   install the first term of mij; the space of fSij in used in the process

    forTerm       = -tField(E11,ijkN)/ffVar(CV_TURB_DENS,ijkN)* &
                     modStrain(ijkN)

    mij(E11,ijkN) = forTerm*fSij(E11,ijkN)
    mij(E12,ijkN) = forTerm*fSij(E12,ijkN)
    mij(E13,ijkN) = forTerm*fSij(E13,ijkN)
    mij(E22,ijkN) = forTerm*fSij(E22,ijkN)
    mij(E23,ijkN) = forTerm*fSij(E23,ijkN)
    mij(E33,ijkN) = forTerm*fSij(E33,ijkN)

! - as final step we determine contribution of test-filtered f-level term;
!   note we have bar[rho] available at cell faces bRho stored in fVar; 
!   in fact we have 1/bRho in fVar(CV_TURB_DENS,:), so we devide by bRho
!   instead of multiply and compute modulus of Sij: |S(u)|

    modStrain(ijkN) = sqrt(sij(E11,ijkN)*sij(E11,ijkN) &
                          +sij(E22,ijkN)*sij(E22,ijkN) &
                          +sij(E11,ijkN)*sij(E22,ijkN) &
                          +sij(E12,ijkN)*sij(E12,ijkN) &
                          +sij(E13,ijkN)*sij(E13,ijkN) &
                          +sij(E23,ijkN)*sij(E23,ijkN))

! - finish computation of mij; filterwidth^2 at f-level is tField/kappa2;
!   build terms need to be filtered

    forTerm = tField(E11,ijkN)/fVar(CV_TURB_DENS,ijkN)*rkappa2* &
              modStrain(ijkN)

    field(E11,ijkN) =  forTerm*sij(E11,ijkN)
    field(E12,ijkN) =  forTerm*sij(E12,ijkN)
    field(E13,ijkN) =  forTerm*sij(E13,ijkN)
    field(E22,ijkN) =  forTerm*sij(E22,ijkN)
    field(E23,ijkN) =  forTerm*sij(E23,ijkN)
    field(E33,ijkN) = -forTerm*(sij(E11,ijkN)+sij(E22,ijkN))
  ENDDO

! perform test filtering:

  idBeg = E11
  idEnd = E33
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,tNDel,idBeg,idEnd,field,tField )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,tNDel,idBeg,idEnd,field,tField )
  ENDIF
#endif
#ifdef RFLU
!  CALL RFLU_InterpFaces2Faces( region,field,tField )
  tField(idBeg:idEnd,:) = field(idBeg:idEnd,:)
#endif

! and complete components of mij:

  DO ijkN=ibn,ien
    mij(E11,ijkN) = mij(E11,ijkN) + tField(E11,ijkN)
    mij(E12,ijkN) = mij(E12,ijkN) + tField(E12,ijkN)
    mij(E13,ijkN) = mij(E13,ijkN) + tField(E13,ijkN)
    mij(E22,ijkN) = mij(E22,ijkN) + tField(E22,ijkN)
    mij(E23,ijkN) = mij(E23,ijkN) + tField(E23,ijkN)
    mij(E33,ijkN) = mij(E33,ijkN) + tField(E33,ijkN)
  ENDDO

! deallocate temporary arrays

  DEALLOCATE( modStrain,field,tField )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesMij

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesMij.F90,v $
! Revision 1.8  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2005/12/31 23:29:02  wasistho
! copy in lack of filter in rocflu temporarily
!
! Revision 1.5  2005/12/28 03:24:39  wasistho
! simply copy variables in waiting for Rocflu's FF-interpolation
!
! Revision 1.4  2004/05/28 02:01:12  wasistho
! update unstructured grid LES
!
! Revision 1.3  2004/03/19 02:50:37  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************








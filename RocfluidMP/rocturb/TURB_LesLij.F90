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
!        ijk    = ijk-face is being treated
!        nDel   = three components filter width paramater
!
! Output: lij = resulting stress tensor of scale similarity or Leonard term
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesLij.F90,v 1.9 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesLij( region,ijk,nDel,lij )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLU
  USE RFLU_ModInterpolation, ONLY : RFLU_InterpCells2Faces
#endif
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetNodeOffset
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
  INTEGER                 :: ijk, nDel(DIRI:DIRK)
  REAL(RFREAL), POINTER   :: lij(:,:)

! ... loop variables
  INTEGER :: ijkN

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ibn,ien,idBeg,idEnd
  REAL(RFREAL), POINTER :: fVar(:,:),ffVar(:,:)
  REAL(RFREAL), POINTER :: tens(:,:),tensBar(:,:)
#ifdef RFLO
  INTEGER :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER :: iLev,iNOff,ijNOff
#endif
#ifdef RFLU
  REAL(RFREAL), POINTER :: cv(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesLij.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesLij',&
  'TURB_LesLij.F90' )

! get indices and pointers --------------------------------------------------

#ifdef RFLO
  iLev   = region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  fVar  => region%levels(iLev)%turb%fVar
  ffVar => region%levels(iLev)%turb%ffVar
  ibn   =  IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien   =  IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
#endif
#ifdef RFLU
  cv    => region%mixt%cv
  fVar  => region%turb%fVar
  ffVar => region%turb%ffVar
  ibn   =  1
  ien   =  region%grid%nFaces
#endif

! allocate temporary arrays

  ALLOCATE( tens(E11:E33,ibn:ien),tensBar(E11:E33,ibn:ien) )

! ff-filtering bRhoF and bRhoUiF (i=1,3) stored in fVar to get 
! bar[bRhoF] and bar[bRhoUiF] stored in fbVar

  idBeg = CV_TURB_DENS
  idEnd = CV_TURB_ZMOM
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,nDel,idBeg,idEnd,fVar,ffVar )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,nDel,idBeg,idEnd,fVar,ffVar )
  ENDIF
#endif
#ifdef RFLU
  CALL RFLU_InterpCells2Faces( region,cv(idBeg:idEnd,:),ffVar )
#endif

! for efficiency we store 1/bRhoF in bRhoF and 1/bbRhoF in bbRhoF;
! next, build components of scale-similarity (Bardina) turbulent stress 
! tensor lij-component: store bRhoUi*bRhoUj/bRho in tens

  DO ijkN = ibn,ien
    fVar(CV_TURB_DENS,ijkN)=1._RFREAL/fVar(CV_TURB_DENS,ijkN)
    ffVar(CV_TURB_DENS,ijkN)=1._RFREAL/ffVar(CV_TURB_DENS,ijkN)

    tens(E11,ijkN) = fVar(CV_TURB_XMOM,ijkN)*fVar(CV_TURB_XMOM,ijkN)* &
                     fVar(CV_TURB_DENS,ijkN)
    tens(E12,ijkN) = fVar(CV_TURB_XMOM,ijkN)*fVar(CV_TURB_YMOM,ijkN)* &
                     fVar(CV_TURB_DENS,ijkN)
    tens(E13,ijkN) = fVar(CV_TURB_XMOM,ijkN)*fVar(CV_TURB_ZMOM,ijkN)* &
                     fVar(CV_TURB_DENS,ijkN)
    tens(E22,ijkN) = fVar(CV_TURB_YMOM,ijkN)*fVar(CV_TURB_YMOM,ijkN)* &
                     fVar(CV_TURB_DENS,ijkN)
    tens(E23,ijkN) = fVar(CV_TURB_YMOM,ijkN)*fVar(CV_TURB_ZMOM,ijkN)* &
                     fVar(CV_TURB_DENS,ijkN)
    tens(E33,ijkN) = fVar(CV_TURB_ZMOM,ijkN)*fVar(CV_TURB_ZMOM,ijkN)* &
                     fVar(CV_TURB_DENS,ijkN)
  ENDDO

! ff-filter bRhoUi*bRhoUj/bRho stored in tens to get 
! bar[bRhoUi*bRhoUj/bRho] in tensBar

  idBeg = E11
  idEnd = E33
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,nDel,idBeg,idEnd,tens,tensBar )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,nDel,idBeg,idEnd,tens,tensBar )
  ENDIF
#endif
#ifdef RFLU
!  CALL RFLU_InterpFaces2Faces( region,tens,tensBar )
  tensBar(idBeg:idEnd,:) = tens(idBeg:idEnd,:)
#endif

! combine and find lij=bar[bRhoUi*bRhoUj/bRho]-bbRhoUi*bbRhoUj/bbRho

  DO ijkN=ibn,ien
    lij(E11,ijkN)=tensBar(E11,ijkN)-ffVar(CV_TURB_DENS,ijkN)* &
                  ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_XMOM,ijkN)

    lij(E12,ijkN)=tensBar(E12,ijkN)-ffVar(CV_TURB_DENS,ijkN)* &
                  ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_YMOM,ijkN)

    lij(E13,ijkN)=tensBar(E13,ijkN)-ffVar(CV_TURB_DENS,ijkN)* &
                  ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)

    lij(E22,ijkN)=tensBar(E22,ijkN)-ffVar(CV_TURB_DENS,ijkN)* &
                  ffVar(CV_TURB_YMOM,ijkN)*ffVar(CV_TURB_YMOM,ijkN)

    lij(E23,ijkN)=tensBar(E23,ijkN)-ffVar(CV_TURB_DENS,ijkN)* &
                  ffVar(CV_TURB_YMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)

    lij(E33,ijkN)=tensBar(E33,ijkN)-ffVar(CV_TURB_DENS,ijkN)* &
                  ffVar(CV_TURB_ZMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)
  ENDDO

! deallocate temporary arrays

  DEALLOCATE( tens,tensBar )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesLij

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesLij.F90,v $
! Revision 1.9  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2005/12/31 23:28:56  wasistho
! copy in lack of filter in rocflu temporarily
!
! Revision 1.6  2005/12/28 03:24:34  wasistho
! simply copy variables in waiting for Rocflu's FF-interpolation
!
! Revision 1.5  2004/08/07 01:06:57  wasistho
! bugfixed, defined pointer cv
!
! Revision 1.4  2004/05/28 02:01:00  wasistho
! update unstructured grid LES
!
! Revision 1.3  2004/03/19 02:50:31  wasistho
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








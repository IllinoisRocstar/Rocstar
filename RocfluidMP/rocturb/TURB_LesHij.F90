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
! Purpose: Get the Leonard stress of the dynamic mixed model.
!
! Description: In this subroutine we determine the H_{ij} components in the
!              dynamic mixed formulation, containing the scale similarity 
!              contribution. We assume L_{ij} to be given already on the cell 
!              face and we generate and subtract Hij from it. On exit, Hij
!              is stored in Lij.
!
! Input: region = data of current region
!        ijk    = ijk-face is being treated
!
! Output: Hij term in dynamic mixed formulation.
!
! Notes: Hij is stored in turb%lij.
!
!******************************************************************************
!
! $Id: TURB_LesHij.F90,v 1.8 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesHij( region,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_LesLij
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesUniFiltFF, TURB_FloLesGenFiltFF, & 
                                 TURB_FloLesGenC2f

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
  INTEGER                 :: ijk

! ... loop variables
  INTEGER :: i, j, k, ijkN

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ibn,ien,idBeg,idEnd
  INTEGER :: NDel(DIRI:DIRK),tNDel(DIRI:DIRK)
  REAL(RFREAL), POINTER :: fVar(:,:),ffVar(:,:),lij(:,:),field(:,:)
  REAL(RFREAL), POINTER :: hij(:,:)
#ifdef RFLO
  INTEGER :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER :: iLev,iNOff,ijNOff
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesHij.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesHij',&
  'TURB_LesHij.F90' )

! get indices and pointers --------------------------------------------------

#ifdef RFLO
  iLev   = region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  fVar  => region%levels(iLev)%turb%fVar
  ffVar => region%levels(iLev)%turb%ffVar
  lij   => region%levels(iLev)%turb%lij

  IF (ijk==DIRI) THEN
    field => region%levels(iLev)%turb%fISij
  ELSEIF (ijk==DIRJ) THEN 
    field => region%levels(iLev)%turb%fJSij
  ELSEIF (ijk==DIRK) THEN 
    field => region%levels(iLev)%turb%fKSij
  ENDIF
  ibn   = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien   = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
#endif

! allocate temporary arrays

  ALLOCATE( hij(E11:E33,ibn:ien) )

! test filter width is twice bar filter width

  nDel(DIRI) = region%turbInput%filterWidth(DIRI)
#ifdef RFLO
  nDel(DIRJ) = region%turbInput%filterWidth(DIRJ)
  nDel(DIRK) = region%turbInput%filterWidth(DIRK)
#endif

  tNDel(DIRI)=2*nDel(DIRI)
#ifdef RFLO
  tNDel(DIRJ)=2*nDel(DIRJ)
  tNDel(DIRK)=2*nDel(DIRK)
#endif

! for efficiency we keep 1/tbRho in tbRho, leftover from lesLij;
! with this information we can generate the first term in Hij
! hij-component; store tbRhoUi*tbRhoUj/tbRho in hij:

  DO ijkN = ibn,ien
    hij(E11,ijkN)=ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_XMOM,ijkN)* &
                  ffVar(CV_TURB_DENS,ijkN)
    hij(E12,ijkN)=ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_YMOM,ijkN)* &
                  ffVar(CV_TURB_DENS,ijkN)
    hij(E13,ijkN)=ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)* &
                  ffVar(CV_TURB_DENS,ijkN)
    hij(E22,ijkN)=ffVar(CV_TURB_YMOM,ijkN)*ffVar(CV_TURB_YMOM,ijkN)* &
                  ffVar(CV_TURB_DENS,ijkN)
    hij(E23,ijkN)=ffVar(CV_TURB_YMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)* &
                  ffVar(CV_TURB_DENS,ijkN)
    hij(E33,ijkN)=ffVar(CV_TURB_ZMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)* &
                  ffVar(CV_TURB_DENS,ijkN)
  ENDDO

! and we double-filter this, first with bar filter

  idBeg = E11
  idEnd = E33
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,nDel,idBeg,idEnd,hij,field )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,nDel,idBeg,idEnd,hij,field )
  ENDIF
#endif

! then with test-filter:

  idBeg = E11
  idEnd = E33
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,tNDel,idBeg,idEnd,field,hij )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,tNDel,idBeg,idEnd,field,hij )
  ENDIF
#endif

! subtract this first term component of hij from lij
! and we put 1/tbRho back to tbRho for the next treatment

  DO ijkN = ibn,ien
    lij(E11,ijkN) = lij(E11,ijkN) - hij(E11,ijkN)
    lij(E12,ijkN) = lij(E12,ijkN) - hij(E12,ijkN)
    lij(E13,ijkN) = lij(E13,ijkN) - hij(E13,ijkN)
    lij(E22,ijkN) = lij(E22,ijkN) - hij(E22,ijkN)
    lij(E23,ijkN) = lij(E23,ijkN) - hij(E23,ijkN)
    lij(E33,ijkN) = lij(E33,ijkN) - hij(E33,ijkN)
    ffVar(CV_TURB_DENS,ijkN) = 1._RFREAL/ffVar(CV_TURB_DENS,ijkN)
  ENDDO

! next, we build second term of the first component in hij for which 
! we apply bar-filter and test-filter to test-bar filtered fields;
! first we filter tbRho and tbRhoUi with bar-filter and store in fVar

  idBeg = CV_TURB_DENS
  idEnd = CV_TURB_ZMOM
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,nDel,idBeg,idEnd,ffVar,fVar )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,nDel,idBeg,idEnd,ffVar,fVar )
  ENDIF
#endif

! then apply test filter and store in ffVar again

#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) then
    CALL TURB_FloLesUniFiltFF( region,ijk,tNDel,idBeg,idEnd,fVar,ffVar )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,tNDel,idBeg,idEnd,fVar,ffVar )
  ENDIF
#endif

! we can now finish the first contribution in hij and subtract it from lij
! for efficiency we store 1/tb(tbRho) in tb(tbRho)

  DO ijkN = ibn,ien
    ffVar(CV_TURB_DENS,ijkN) = 1._RFREAL/ffVar(CV_TURB_DENS,ijkN)
    lij(E11,ijkN) = lij(E11,ijkN) + ffVar(CV_TURB_DENS,ijkN)* &
                    ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_XMOM,ijkN)

    lij(E12,ijkN) = lij(E12,ijkN) + ffVar(CV_TURB_DENS,ijkN)* &
                    ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_YMOM,ijkN)

    lij(E13,ijkN) = lij(E13,ijkN) + ffVar(CV_TURB_DENS,ijkN)* &
                    ffVar(CV_TURB_XMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)

    lij(E22,ijkN) = lij(E22,ijkN) + ffVar(CV_TURB_DENS,ijkN)* &
                    ffVar(CV_TURB_YMOM,ijkN)*ffVar(CV_TURB_YMOM,ijkN)

    lij(E23,ijkN) = lij(E23,ijkN) + ffVar(CV_TURB_DENS,ijkN)* &
                    ffVar(CV_TURB_YMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)

    lij(E33,ijkN) = lij(E33,ijkN) + ffVar(CV_TURB_DENS,ijkN)* &
                    ffVar(CV_TURB_ZMOM,ijkN)*ffVar(CV_TURB_ZMOM,ijkN)
  ENDDO

! we now come to the second component which involves the test-filtered
! Bardina term; first recompute density and velocities at faces and 
! generate the Bardina term and put the result in field

#ifdef RFLO
  CALL TURB_FloLesGenC2F( region,ijk )
#endif

  CALL TURB_LesLij( region,ijk,nDel,field )

! test-filter Bardina component stored in field and put in hij

  idBeg = E11
  idEnd = E33
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) THEN
    CALL TURB_FloLesUniFiltFF( region,ijk,tNDel,idBeg,idEnd,field,hij )
  ELSE
    CALL TURB_FloLesGenFiltFF( region,ijk,tNDel,idBeg,idEnd,field,hij )
  ENDIF
#endif

! complete hij and subtract it from lij

  DO ijkN = ibn,ien
    lij(E11,ijkN) = lij(E11,ijkN) + hij(E11,ijkN)
    lij(E12,ijkN) = lij(E12,ijkN) + hij(E12,ijkN)
    lij(E13,ijkN) = lij(E13,ijkN) + hij(E13,ijkN)
    lij(E22,ijkN) = lij(E22,ijkN) + hij(E22,ijkN)
    lij(E23,ijkN) = lij(E23,ijkN) + hij(E23,ijkN)
    lij(E33,ijkN) = lij(E33,ijkN) + hij(E33,ijkN)
  ENDDO

! deallocate temporary arrays

  DEALLOCATE( hij )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesHij

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesHij.F90,v $
! Revision 1.8  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2004/07/30 22:35:11  wasistho
! replaced floLesUniC2F by floLesGenC2F
!
! Revision 1.5  2004/05/28 02:01:06  wasistho
! update unstructured grid LES
!
! Revision 1.4  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.3  2004/03/19 02:50:25  wasistho
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








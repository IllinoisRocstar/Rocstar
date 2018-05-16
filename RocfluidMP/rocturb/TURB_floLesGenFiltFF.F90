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
! Purpose: Perform face to face filtering.
!
! Description: The filtering is performed in i, j and k direction, subsequently.
!              For each direction, distinction is made between homogeneous
!              and non-homogeneous direction. If the filtering direction is
!              homogeneous, uniform filtering is carried out, otherwise
!              nonuniform filtering.
!
! Input: region  = data of current region 
!        ijk     = ijk-face being treated
!        nDel    = three components filter width parameter
!        idBeg   = begin variable index to be filtered
!        idEnd   = end variable index to be filtered
!        fVar    = face variable to be filtered
!
! Output: fbVar  = filtered face variable
!
! Notes: This routine is only relevant if non-uniform filter is selected.
!
!******************************************************************************
!
! $Id: TURB_floLesGenFiltFF.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenFiltFF( region,ijk,nDel,idBeg,idEnd,fVar,fbVar )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesUniFiltFFI, &
                               TURB_FloLesUniFiltFFJ, TURB_FloLesUniFiltFFK, &
                               TURB_FloLesGenFiltFFI, TURB_FloLesGenFiltFFJ, &
                               TURB_FloLesGenFiltFFK, TURB_FloExtrapolFaceVec
  USE ModTurbulence
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region)        :: region
  INTEGER               :: ijk,nDel(DIRI:DIRK),idBeg,idEnd
  REAL(RFREAL), POINTER :: fVar(:,:),fbVar(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, ijkN

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER           :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER           :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER           :: iLev,iNOff,ijNOff,ibn,ien
  INTEGER           :: homDir(DIRI:DIRK)

  REAL(RFREAL)      :: fact1(FILWIDTH_FOUR), fact2(FILWIDTH_FOUR)
  REAL(RFREAL), POINTER :: wrkbar(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenFiltFF.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenFiltFF',&
  'TURB_floLesGenFiltFF.F90' )

! get indices, parameters and pointers ----------------------------------------

  iLev      = region%currLevel
  homDir(:) = region%turbInput%homDir(:)
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  ALLOCATE( wrkbar(idBeg:idEnd,ibn:ien) )

  fact1(FILWIDTH_ONE)  = 0.125_RFREAL
  fact2(FILWIDTH_ONE)  = 0.75_RFREAL
  fact1(FILWIDTH_TWO)  = 0.25_RFREAL
  fact2(FILWIDTH_TWO)  = 0.50_RFREAL
  fact1(FILWIDTH_FOUR) = 0.125_RFREAL
  fact2(FILWIDTH_FOUR) = 0.25_RFREAL

!------------------------------------------------------------------------
! we begin with integration over I-direction

  ibeg = idnbeg+1
  iend = idnend-1
  jbeg = jdnbeg
  jend = jdnend
  kbeg = kdnbeg
  kend = kdnend

  IF (homDir(DIRI) /= OFF) THEN
    CALL TURB_FloLesUniFiltFFI( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                             ijNOff,nDel,idBeg,idEnd,fact1,fact2,fVar,fbVar )
  ELSE
    CALL TURB_FloLesGenFiltFFI( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                             iNOff,ijNOff,nDel,idBeg,idEnd,fVar,fbVar )
  ENDIF

! extrapolate filtered variables to outers ijk-face layers

  CALL TURB_FloExtrapolFaceVec( region,DIRI,idBeg,idEnd,fbVar )

!------------------------------------------------------------------
! next, integrate over J-direction

  ibeg = idnbeg
  iend = idnend
  jbeg = jdnbeg+1
  jend = jdnend-1
  kbeg = kdnbeg
  kend = kdnend

  IF (homDir(DIRJ) /= OFF) THEN
    CALL TURB_FloLesUniFiltFFJ( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                             ijNOff,nDel,idBeg,idEnd,fact1,fact2,fbVar,wrkBar )
  ELSE
    CALL TURB_FloLesGenFiltFFJ( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                             iNOff,ijNOff,nDel,idBeg,idEnd,fbVar,wrkBar )
  ENDIF

! extrapolate filtered variables to outerst ijk-face layers

  CALL TURB_FloExtrapolFaceVec( region,DIRJ,idBeg,idEnd,wrkBar )

!------------------------------------------------------------------
! finally, integrate over K-direction

  ibeg = idnbeg
  iend = idnend
  jbeg = jdnbeg
  jend = jdnend
  kbeg = kdnbeg+1
  kend = kdnend-1

  IF (homDir(DIRK) /= OFF) THEN
    CALL TURB_FloLesUniFiltFFK( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                             ijNOff,nDel,idBeg,idEnd,fact1,fact2,wrkBar,fbVar )
  ELSE
    CALL TURB_FloLesGenFiltFFK( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                             iNOff,ijNOff,nDel,idBeg,idEnd,wrkBar,fbVar )
  ENDIF

! extrapolate filtered variables to outerst ijk-face layers

  CALL TURB_FloExtrapolFaceVec( region,DIRK,idBeg,idEnd,fbVar )

! deallocate temporary arrays

  DEALLOCATE( wrkBar )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenFiltFF

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenFiltFF.F90,v $
! Revision 1.4  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************








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
! Purpose: contracting tensors M_{ij} and L_{ij}
!
! Description: We determine the dynamic coefficient by contracting tensors 
!              M_{ij} and (L_{ij}-H_{ij}) and evalutating 
!              C=<M_{ij}(L_{ij}-H_{ij})>/<M_{ij}M_{ij}> in which <.> denotes 
!              an averaging over the homogeneous directions. 
!              On entry lij contains L_{ij}-H_{ij} with H_{ij} appropriate 
!              for the present model. A 1D vector representing model coef is 
!              returned and it is assumed that all information is available 
!              at cell faces.
!
! Input: region  = data of current region
!        ijk     = ijk-face is being treated
!
! Output: coef = turb%coef = model coefficient
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesContract.F90,v 1.9 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesContract( region,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLU
  USE ModBndPatch, ONLY   : t_patch
#endif
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesAverageFace

#include "Indexing.h"
#endif
  USE TURB_ModInterfaces, ONLY : TURB_StatFCollector
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
  INTEGER :: iN, iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_patch), POINTER  :: patch
#endif

  INTEGER :: ibn, ien, nDel(DIRI:DIRK)
  REAL(RFREAL), POINTER :: coef(:,:),mij(:,:),lij(:,:), colVar(:,:)
  REAL(RFREAL), POINTER :: mijMij(:),mijLij(:),avgMMij(:),avgMLij(:)

#ifdef RFLO
  INTEGER :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER :: iLev,iNOff,ijNOff
#endif
#ifdef RFLU
  INTEGER :: nPatches, nBFaces
  REAL(RFREAL), POINTER :: bCoef(:,:),bMij(:,:),bLij(:,:),colBVar(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesContract.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesContract',&
  'TURB_LesContract.F90' )

! get indices and pointers --------------------------------------------------

  nDel(:) =  region%turbInput%filterWidth(:)

#ifdef RFLO
  iLev    =  region%currLevel

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  coef    => region%levels(ilev)%turb%coef
  lij     => region%levels(ilev)%turb%lij
  mij     => region%levels(ilev)%turb%mij
#endif
#ifdef RFLU
  ibn = 1
  ien = region%grid%nFaces

  coef    => region%turb%coef
  lij     => region%turb%lij
  mij     => region%turb%mij
#endif

! allocate temporary arrays

  ALLOCATE( mijMij(ibn:ien) )
  ALLOCATE( mijLij(ibn:ien) )
  ALLOCATE( avgMMij(ibn:ien) )
  ALLOCATE( avgMLij(ibn:ien) )

! first contract the numerator and denominator

  DO iN = ibn,ien

! - M_{ij}*M_{ij}:

    mijMij(iN) = 2._RFREAL*(mij(E11,iN)*mij(E11,iN) &
                           +mij(E22,iN)*mij(E22,iN) &
                           +mij(E11,iN)*mij(E22,iN) &
                           +mij(E12,iN)*mij(E12,iN) &
                           +mij(E13,iN)*mij(E13,iN) &
                           +mij(E23,iN)*mij(E23,iN))

! - M_{ij}*(L_{ij}-H_{ij}):

    mijLij(iN) = mij(E11,iN)*lij(E11,iN)  &
                +mij(E22,iN)*lij(E22,iN)  &
                +mij(E33,iN)*lij(E33,iN)  &
          +2.0d0*mij(E12,iN)*lij(E12,iN)  &
          +2.0d0*mij(E13,iN)*lij(E13,iN)  &
          +2.0d0*mij(E23,iN)*lij(E23,iN)
  ENDDO

! dynamic coefficient is obtained by averaging numerator and denominator
! first, denominator

#ifdef RFLO
  CALL TURB_FloLesAverageFace( region,ijk,mijMij,avgMMij )
#endif
#ifdef RFLU
  avgMMij = mijMij
#endif

! then, numerator:

#ifdef RFLO
  CALL TURB_FloLesAverageFace( region,ijk,mijLij,avgMLij )
#endif
#ifdef RFLU
  avgMLij = mijLij
#endif

! we assign avgMLij/avgMMij to coef in which we also limit from below with 0; 
! moreover, if all the filterwidth is zero coef set to zero

  IF (nDel(DIRI)==0 .AND. nDel(DIRJ)==0 .AND. nDel(DIRK)==0) THEN
    coef=0.0_RFREAL
  ELSE
    DO iN = ibn,ien
      IF (avgMLij(iN) > 0.1_RFREAL*avgMMij(iN)) avgMLij(iN) = 0._RFREAL
      coef(1,iN)= MAX(0._RFREAL,avgMLij(iN)/(avgMMij(iN)+REAL_SMALL))
    ENDDO
  ENDIF

!#ifdef STATS
! if desired, collect quantities of interest at cell centers for statistics

!  IF (region%turbInput%nSt > 0) THEN
!    ALLOCATE( colVar(2,ibn:ien) )
!    colVar(1,:) = avgMMij
!    colVar(2,:) = avgMLij
!#ifdef RFLO
!    CALL TURB_StatFCollector( region,ijk,1,2,colVar )
!    DEALLOCATE( colVar )
!#endif
!  ENDIF
!#endif

! deallocate retired arrays
  DEALLOCATE( mijMij,mijLij,avgMMij,avgMLij )

! RFLU boundary treatments =====================================================
#ifdef RFLU

  nPatches = region%grid%nPatches
  nBFaces  = 0

  DO iPatch = 1,nPatches
    patch   => region%patches(iPatch)
    nBFaces =  nBFaces + patch%nBTris + patch%nBQuads
  END DO ! iPatch
  ibn = 1
  ien = nBFaces

  bCoef   => region%turb%bCoef
  bLij    => region%turb%bLij
  bMij    => region%turb%bMij

! allocate temporary arrays

  ALLOCATE( mijMij(ibn:ien) )
  ALLOCATE( mijLij(ibn:ien) )
  ALLOCATE( avgMMij(ibn:ien) )
  ALLOCATE( avgMLij(ibn:ien) )

! first contract the numerator and denominator

  DO iN = ibn,ien

! - M_{ij}*M_{ij}:

    mijMij(iN) = 2._RFREAL*(bMij(E11,iN)*bMij(E11,iN) &
                           +bMij(E22,iN)*bMij(E22,iN) &
                           +bMij(E11,iN)*bMij(E22,iN) &
                           +bMij(E12,iN)*bMij(E12,iN) &
                           +bMij(E13,iN)*bMij(E13,iN) &
                           +bMij(E23,iN)*bMij(E23,iN))

! - M_{ij}*(L_{ij}-H_{ij}):

    mijLij(iN) = bMij(E11,iN)*bLij(E11,iN)  &
                +bMij(E22,iN)*bLij(E22,iN)  &
                +bMij(E33,iN)*bLij(E33,iN)  &
          +2.0d0*bMij(E12,iN)*bLij(E12,iN)  &
          +2.0d0*bMij(E13,iN)*bLij(E13,iN)  &
          +2.0d0*bMij(E23,iN)*bLij(E23,iN)
  ENDDO

! dynamic coefficient is obtained by averaging numerator and denominator
! first, denominator

  avgMMij = mijMij

! then, numerator:

  avgMLij = mijLij

! we assign avgMLij/avgMMij to coef in which we also limit from below with 0; 
! moreover, if all the filterwidth is zero coef set to zero

  IF (nDel(DIRI)==0 .AND. nDel(DIRJ)==0 .AND. nDel(DIRK)==0) THEN
    bCoef=0.0_RFREAL
  ELSE
    DO iN = ibn,ien
      IF (avgMLij(iN) > 0.1_RFREAL*avgMMij(iN)) avgMLij(iN) = 0._RFREAL
      bCoef(1,iN)= MAX(0._RFREAL,avgMLij(iN)/(avgMMij(iN)+REAL_SMALL))
    ENDDO
  ENDIF

!#ifdef STATS
! if desired, collect quantities of interest at cell centers for statistics

!  IF (region%turbInput%nSt > 0) THEN
!    ALLOCATE( colBVar(2,ibn:ien) )
!    colBVar(1,:) = avgMMij
!    colBVar(2,:) = avgMLij
!    CALL TURB_StatFCollector( region,ijk,1,2,colVar,colBVar )
!    DEALLOCATE( colVar,colBVar )
!  ENDIF
!#endif

! deallocate retired arrays
  DEALLOCATE( mijMij,mijLij,avgMMij,avgMLij )

#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesContract

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesContract.F90,v $
! Revision 1.9  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/01/30 23:12:08  wasistho
! removed comment after endif RFLU
!
! Revision 1.6  2004/10/22 23:18:26  wasistho
! comment out statistics collection of Mij and Lij
!
! Revision 1.5  2004/05/28 02:01:35  wasistho
! update unstructured grid LES
!
! Revision 1.4  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.3  2004/03/19 02:49:35  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.4  2004/02/18 21:30:26  wasistho
! set <Mij.Lij> to zero when Lij/Mij is excessively large
!
! Revision 1.3  2003/05/24 02:07:33  wasistho
! turbulence statistics expanded
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************








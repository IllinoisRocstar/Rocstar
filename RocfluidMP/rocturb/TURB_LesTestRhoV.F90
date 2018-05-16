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
! Purpose: Compute test filtered density and velocities at cell centers.
!
! Description: We filter mixture cv using cell to cell (test) filtering.
!              The cell filtered cv is then processed to get Favre filtered
!              velocities.
!
! Input: region = data of current region
!
! Output: ccVar = cell filtered rho and u_i
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesTestRhoV.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesTestRhoV( region,ibc,iec )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE TURB_ModInterfaces, ONLY : TURB_FloLesUniFiltCC, TURB_FloLesGenFiltCC 

#include "Indexing.h"
#endif
  USE ModTurbulence
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: ibc, iec

! ... loop variables
  INTEGER :: i, j, k, ijkC

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: idBeg,idEnd,tNDel(DIRI:DIRK)
  REAL(RFREAL), POINTER :: cv(:,:), ccVar(:,:)

#ifdef RFLO
  INTEGER :: iLev
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesTestRhoV.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesTestRhoV',&
  'TURB_LesTestRhoV.F90' )

! get indices and pointers --------------------------------------------------

#ifdef RFLO
  iLev  =  region%currLevel
  cv    => region%levels(iLev)%mixt%cv
  ccVar => region%levels(iLev)%turb%ccVar
#endif
#ifdef RFLU
  cv    => region%mixt%cv
  ccVar => region%turb%ccVar
#endif

! test filter width is twice bar filter width

  tNDel(DIRI) = 2*region%turbInput%filterWidth(DIRI)
#ifdef RFLO
  tNDel(DIRJ) = 2*region%turbInput%filterWidth(DIRJ)
  tNDel(DIRK) = 2*region%turbInput%filterWidth(DIRK)
#endif

! we calculate vi=test[bar(rho*ui)]/test[bar(rho)] from bar(rho) and 
! bar(rho*ui) already available as mixt%cv; the filtered conservative
! variables and returned at cell centers including dummies in ccVar

  idBeg = CV_TURB_DENS
  idEnd = CV_TURB_ZMOM
#ifdef RFLO
  IF (region%turbInput%filterType == FILTYPE_UNIFORM) THEN
    CALL TURB_FloLesUniFiltCC( region,tNDel,idBeg,idEnd,cv,ccVar )
  ELSE
    CALL TURB_FloLesGenFiltCC( region,tNDel,idBeg,idEnd,cv,ccVar )
  ENDIF
#endif
#ifdef RFLU
!  CALL TURB_FluLesFiltCC( region,tNDel,idBeg,idEnd,cv,ccVar )
  ccVar(idBeg:idEnd,:) = cv(idBeg:idEnd,:)
#endif

! for efficiency we store 1/test(bar[rho]) in ccVar(CV_TURB_DENS,:)

  DO ijkC=ibc,iec
    ccVar(CV_TURB_DENS,ijkC) = 1._RFREAL/ccVar(CV_TURB_DENS,ijkC)
  ENDDO

! then compute vi and store these in ccVar

  DO ijkC=ibc,iec
    ccVar(CV_TURB_UVEL,ijkC) = ccVar(CV_TURB_XMOM,ijkC)* &
                               ccVar(CV_TURB_DENS,ijkC)
    ccVar(CV_TURB_VVEL,ijkC) = ccVar(CV_TURB_YMOM,ijkC)* &
                               ccVar(CV_TURB_DENS,ijkC)
    ccVar(CV_TURB_WVEL,ijkC) = ccVar(CV_TURB_ZMOM,ijkC)* &
                               ccVar(CV_TURB_DENS,ijkC)
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesTestRhoV

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesTestRhoV.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/12/31 23:28:15  wasistho
! copy cv to ccVar for rocflu temporarily
!
! Revision 1.3  2004/03/19 02:51:34  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.3  2003/05/16 05:43:44  wasistho
! modified array range of CC-filtered
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************








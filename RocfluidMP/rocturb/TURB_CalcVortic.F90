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
! Purpose: Compute vorticity components at cells.
!
! Description: Velocity gradients at faces first averaged into cells.
!              The resulting cell gradients are then used to compute vorticity
!              components.
!
! Input: region  = data of current region
!
! Output: x,y,z-vorticity, stored in turb%vort(XCOORD:ZCOORD,:)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_CalcVortic.F90,v 1.7 2009/08/26 12:28:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CalcVortic( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ic, idGr

! ... local variables
  TYPE(t_global), POINTER :: global
  
  INTEGER :: errorFlag
  REAL(RFREAL), POINTER :: vort(:,:)

#ifdef RFLO
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: ijkC, ijkN, ijkNi, ijkNj, ijkNk
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  REAL(RFREAL) :: oo6, grVel(GR_MIXT_UX:GR_MIXT_TZ)
  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
#endif
#ifdef RFLU
  REAL(RFREAL) :: grVel(XCOORD:ZCOORD, GRC_MIXT_XVEL:GRC_MIXT_ZVEL)
  REAL(RFREAL), POINTER :: gradc(:,:,:)
#endif

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'TURB_CalcVortic',&
  'TURB_CalcVortic.F90' )

#ifdef RFLO
! get RFLO dimensions, pointers and parameters --------------------------------
  iLev  =  region%currLevel

  CALL RFLO_GetDimensDummy( region,ilev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  gradi => region%levels(iLev)%mixt%gradi
  gradj => region%levels(iLev)%mixt%gradj
  gradk => region%levels(iLev)%mixt%gradk
  vort  => region%levels(iLev)%turb%vort

  oo6 = 1._RFREAL/6._RFREAL

! compute vorticities

  DO k = kdcbeg,kdcend
    DO j = jdcbeg,jdcend
      DO i = idcbeg,idcend
        ijkC = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        ijkNi= IndIJK(i+1   ,j     ,k     ,iNOff,ijNOff)
        ijkNj= IndIJK(i     ,j+1   ,k     ,iNOff,ijNOff)
        ijkNk= IndIJK(i     ,j     ,k+1   ,iNOff,ijNOff)

        DO idGr = GR_MIXT_UX, GR_MIXT_TZ
          grVel(idGr) = gradi(idGr,ijkN) + gradi(idGr,ijkNi) + &
                        gradj(idGr,ijkN) + gradj(idGr,ijkNj) + &
                        gradk(idGr,ijkN) + gradk(idGr,ijkNk) 
        ENDDO

        vort(XCOORD,ijkC) = oo6*(grVel(GR_MIXT_WY) - grVel(GR_MIXT_VZ))
        vort(YCOORD,ijkC) = oo6*(grVel(GR_MIXT_UZ) - grVel(GR_MIXT_WX))
        vort(ZCOORD,ijkC) = oo6*(grVel(GR_MIXT_VX) - grVel(GR_MIXT_UY))
      ENDDO ! i
    ENDDO   ! j
  ENDDO     ! k
#endif

#ifdef RFLU
! get RFLU dimensions, pointers and parameters --------------------------------

  IF (ASSOCIATED( region%mixt%gradCell ) .eqv. .true.) THEN
    gradc => region%mixt%gradCell
    vort  => region%turb%vort

! - compute vorticities

    DO ic = 1,region%grid%nCellsTot
      vort(XCOORD,ic) = gradc(YCOORD,GRC_MIXT_ZVEL,ic) - &
                        gradc(ZCOORD,GRC_MIXT_YVEL,ic)
      vort(YCOORD,ic) = gradc(ZCOORD,GRC_MIXT_XVEL,ic) - &
                        gradc(XCOORD,GRC_MIXT_ZVEL,ic)
      vort(ZCOORD,ic) = gradc(XCOORD,GRC_MIXT_YVEL,ic) - &
                        gradc(YCOORD,GRC_MIXT_XVEL,ic)
    ENDDO ! ic
  ENDIF   ! gradCell allocated?
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CalcVortic

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_CalcVortic.F90,v $
! Revision 1.7  2009/08/26 12:28:52  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.6  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/16 01:04:49  wasistho
! compute RFLU vorticity only if gradCell allocated
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/19 02:45:02  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2004/02/12 22:16:55  wasistho
! changed RFLO_GetDimensPhys to RFLO_GetDimensDummy
!
! Revision 1.1  2003/08/06 15:58:41  wasistho
! added vorticities computation
!
!
!
!******************************************************************************








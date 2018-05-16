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
! Purpose: update of gradients at the current patch if physical 
!          boundary conditions apply
!
! Description: gradients at physical boundaries computed using half interior 
!              side of control volume; resulting gradients computation from 
!              RFLO_calcGradFaces corrected by taking face averaging
!              at physical boundary for linearly extrapolated variables or
!              multiplicating by two for constantly extrapolated variables
!
! Input: region        = data of current region.
!        patch         = current patch.
!        iBegV, iEndV  = begin and end var index
!        iBegG, iEndG  = begin and end gradient index
!        var           = variables, the gradient of which to be determined.
!  
! Output: gradi, gradj, gradk = gradients at the patch faces subject to 
!                               physical boundary conditions
!
! Notes: two other subroutines contained in this file:
!          RFLO_FinishGradPatchLinear for linearly extrapolated dummy variables
!          RFLO_FinishGradPatchConstant for constantly extrapolated dummy
!
!******************************************************************************
!
! $Id: RFLO_CalcGradPhysBc.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradPhysBc( region,patch,iBegV,iEndV,iBegG,iEndG,var, &
                                gradi,gradj,gradk )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices,RFLO_GetPatchDirection, &
        RFLO_GetCellOffset,RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: iBegV, iEndV, iBegG, iEndG
  REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:), gradk(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, lx, ly, lz

! ... local variables
  INTEGER :: iLev, lbound, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: inode, jnode, knode, idir, jdir, kdir
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkN, ijkN1, ijkC, ijkC1
  INTEGER :: nvar
  INTEGER :: iExtrap_SlipW,iExtrap_Inject

  REAL(RFREAL)          :: sgn, fnx, fny, fnz, avol, rvol, avar(iBegV:iEndV)
  REAL(RFREAL)          :: snx, sny, snz, rsurf
  REAL(RFREAL)          :: sdotgrad(iBegV:iEndV)
  REAL(RFREAL), POINTER :: sFace(:,:), grad(:,:), vol(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradPhysBc',&
  'RFLO_CalcGradPhysBc.F90' )

! get dimensions and pointers-------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound
  bcType = patch%bcType

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nvar = iEndV - iBegV + 1

  vol => region%levels(ilev)%grid%vol

! take the right face vector and make it point outwards ----------------------

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF( lbound==2 .OR. lbound==4 .OR. lbound==6 ) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

  IF (lbound==1 .OR. lbound==2) THEN
    sFace => region%levels(iLev)%grid%si
    grad  => gradi
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    sFace => region%levels(iLev)%grid%sj
    grad  => gradj
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    sFace => region%levels(iLev)%grid%sk
    grad  => gradk
  ENDIF

  IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) &
     iExtrap_SlipW  = patch%mixt%switches(BCSWI_SLIPW_EXTRAP)
  IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) &
     iExtrap_Inject = patch%mixt%switches(BCSWI_INJECT_EXTRAP)

  IF ((bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE &
                           .AND. iExtrap_Slipw==EXTRAPOL_LINEAR) .OR. &
      (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
      (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE &
                            .AND. iExtrap_Inject==EXTRAPOL_LINEAR) .OR. &
      (bcType>=BC_INFLOW .AND. bcType<=BC_INFLOW+BC_RANGE) .OR. &
      (bcType>=BC_OUTFLOW .AND. bcType<=BC_OUTFLOW+BC_RANGE) .OR. &
      (bcType>=BC_FARFIELD .AND. bcType<=BC_FARFIELD+BC_RANGE)) THEN

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          CALL RFLO_FinishGradPatchLinear
        ENDDO
      ENDDO
    ENDDO

  ELSE IF((bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE &
                               .AND. iExtrap_Slipw==EXTRAPOL_CONST ) .OR. &
          (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE &
                                .AND. iExtrap_Inject==EXTRAPOL_CONST )) THEN

    CALL RFLO_FinishGradPatchConstant

  ELSE IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN

! remove the portion normal to the bndry patch --------------------------------
! surface vectors pointed inwards/negative
! 2D loop over patch faces

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          CALL RFLO_FinishGradPatchLinear

          ijkN  = IndIJK(i+inode, j+jnode, k+knode, iNOff,ijNOff)  ! bnd nodes
          rsurf = sFace(XCOORD,ijkN)**2+sFace(YCOORD,ijkN)**2+sFace(ZCOORD,ijkN)**2
          rsurf = 1.0_RFREAL/rsurf
          snx   = sFace(XCOORD,ijkN)*rsurf
          sny   = sFace(YCOORD,ijkN)*rsurf
          snz   = sFace(ZCOORD,ijkN)*rsurf

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar

            sdotgrad(l) = sFace(XCOORD,ijkN)*grad(lx,ijkN)+ &
                          sFace(YCOORD,ijkN)*grad(ly,ijkN)+ &
                          sFace(ZCOORD,ijkN)*grad(lz,ijkN)

            grad(lx,ijkN) = grad(lx,ijkN)-sdotgrad(l)*snx
            grad(ly,ijkN) = grad(ly,ijkN)-sdotgrad(l)*sny
            grad(lz,ijkN) = grad(lz,ijkN)-sdotgrad(l)*snz
          ENDDO

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! everything else 

  ELSE

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

! ==============================================================================
! two patch gradient completion subroutines for linear and constant extrap. case
! ==============================================================================

CONTAINS

  SUBROUTINE RFLO_FinishGradPatchLinear

! correct gradients at boundary whose dummy variables are linearly extrapolated

      ijkC  = IndIJK(i        ,j        ,k        ,iCOff,ijCOff)
      ijkC1 = IndIJK(i-idir   ,j-jdir   ,k-kdir   ,iCOff,ijCOff)
      ijkN  = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)
      ijkN1 = IndIJK(i+inode-idir,j+jnode-jdir,k+knode-kdir,iNOff,ijNOff)

      fnx =  .5_RFREAL*(sFace(XCOORD,ijkN)+sFace(XCOORD,ijkN1))
      fny =  .5_RFREAL*(sFace(YCOORD,ijkN)+sFace(YCOORD,ijkN1))
      fnz =  .5_RFREAL*(sFace(ZCOORD,ijkN)+sFace(ZCOORD,ijkN1))
      avol=  .5_RFREAL*(vol(ijkC)+vol(ijkC1))
      rvol= 2.0_RFREAL/avol  ! - for the devision by half control volume    

      DO l=iBegV,iEndV
        lx=l+iBegG-iBegV 
        ly=lx+nvar
        lz=ly+nvar

! - cancel existing patch contribution for front patches (lbound=1,3,5)

        grad(lx,ijkN)=grad(lx,ijkN)*avol-sgn*fnx*var(l,ijkC1)
        grad(ly,ijkN)=grad(ly,ijkN)*avol-sgn*fny*var(l,ijkC1)
        grad(lz,ijkN)=grad(lz,ijkN)*avol-sgn*fnz*var(l,ijkC1)

! - set correct patch contribution 
!   avg.coeffs. are not explicitly used as they are always 0.5 at linear patch

        avar(l)=.5_RFREAL*(var(l,ijkC)+var(l,ijkC1))
        grad(lx,ijkN)=grad(lx,ijkN)+sgn*fnx*avar(l)
        grad(ly,ijkN)=grad(ly,ijkN)+sgn*fny*avar(l)
        grad(lz,ijkN)=grad(lz,ijkN)+sgn*fnz*avar(l)

! - division by half control volume

        grad(lx,ijkN)=grad(lx,ijkN)*rvol
        grad(ly,ijkN)=grad(ly,ijkN)*rvol
        grad(lz,ijkN)=grad(lz,ijkN)*rvol
      ENDDO

  END SUBROUTINE RFLO_FinishGradPatchLinear

! ##############################################################################

  SUBROUTINE RFLO_FinishGradPatchConstant

! correct gradients at boundary whose dummy variables are constantly extrapolated
! account for half control volume by multiplication by 2

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkN  = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar

            grad(lx,ijkN) = 2.0_RFREAL*grad(lx,ijkN)
            grad(ly,ijkN) = 2.0_RFREAL*grad(ly,ijkN)
            grad(lz,ijkN) = 2.0_RFREAL*grad(lz,ijkN)
          ENDDO

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE RFLO_FinishGradPatchConstant

END SUBROUTINE RFLO_CalcGradPhysBc

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcGradPhysBc.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:37:57  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.9  2004/08/03 22:48:50  wasistho
! changed cell2edge averaging to grid dependent avg
!
! Revision 1.8  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.4  2002/10/15 15:37:59  wasistho
! slipw and injection extrapolation switch
!
! Revision 1.3  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/09/02 22:58:54  wasistho
! RFLO grad routines migrated from rocflo to libflo
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:44:24  wasistho
! made compliant with CODING RULE
!
!******************************************************************************








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
! Purpose: specify gradients at the current patch subject to connecting 
!          boundary condition (block interface or periodic)
!
! Description: if a region side has a patch with connecting bc, treatment
!              applied on whole side to omit complications at inter-patch 
!              junctions; non-conneting patches will be overruled by other bc
!
! Input: region        = data of current region
!        patch         = current patch
!        iConBc        = connecting Bc treatment flag (1:treated, 0: not yet)
!        iBegV, iEndV  = begin and end var index
!        iBegG, iEndG  = begin and end gradient index
!        var           = variables, the gradient of which to be determined
!  
! Output: gradi, gradj, gradk = gradients at the patch faces subject to 
!                connecting boundary conditions
!
! Notes: 1) gradients at edges and corners can be computed more accurately
!           if the required information of face vectors are available
!        2) three other routines contained in this file:
!           - RFLO_SetIndexRange to set index-range to whole region-side
!           - RFLO_CalcSideGrad to compute gradients on region-side
!           - RFLO_CopyPatchEdgeGrad to copy gradients to patch edges
!
!******************************************************************************
!
! $Id: RFLO_CalcGradConnBc.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradConnBc( region,patch,iConBc,iBegV,iEndV,iBegG,iEndG, &
                                var,gradi,gradj,gradk )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: iConBc(6),iBegV,iEndV,iBegG,iEndG
  REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:), gradk(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, lx, ly, lz

! ... local variables
  INTEGER :: iLev, lbound, bcType
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: inbeg, inend, jnbeg, jnend, knbeg, knend
  INTEGER :: idir, jdir, kdir
  INTEGER :: iCOff,ijCOff,iNOff,ijNOff
  INTEGER :: ijkC, im1jkC, ijm1kC, ijkm1C, im1jm1kC, im1jkm1C, ijm1km1C, &
             im1jm1km1C, ijkN, im1jkN, ijm1kN, ijkm1N, nvar

  REAL(RFREAL), POINTER :: aci(:,:), acj(:,:), ack(:,:)
  REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:), vol(:)
  REAL(RFREAL)          :: fnx, fny, fnz, rvol, fface(iBegV:iEndV), tvar

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradConnBc',&
       'RFLO_CalcGradConnBc.F90' )
 
! get cell indices and geometrical parameters --------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound
  bcType = patch%bcType

  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! complete the outer layers --------------------------------------------------

  IF ((bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI+  BC_RANGE) .OR. &
      (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI+  BC_RANGE) .OR. &
      (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE)) THEN

    CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
    CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    nvar = iEndV - iBegV + 1

! - get dimensions and pointers

    aci => region%levels(ilev)%grid%c2eCoI
    acj => region%levels(ilev)%grid%c2eCoJ
    ack => region%levels(ilev)%grid%c2eCoK
    si  => region%levels(iLev)%grid%si
    sj  => region%levels(iLev)%grid%sj
    sk  => region%levels(iLev)%grid%sk
    vol => region%levels(ilev)%grid%vol

    IF (iConBc(lbound)==0) THEN
      IF (lbound==1 .OR. lbound==2) THEN
        jbeg=jpcbeg
        jend=jpcend
        kbeg=kpcbeg
        kend=kpcend
      ELSE IF (lbound==3 .OR. lbound==4) THEN
        ibeg=ipcbeg
        iend=ipcend
        kbeg=kpcbeg
        kend=kpcend
      ELSE IF (lbound==5 .OR. lbound==6) THEN
        ibeg=ipcbeg
        iend=ipcend
        jbeg=jpcbeg
        jend=jpcend
      ENDIF
      CALL RFLO_SetIndexRange
      CALL RFLO_CalcSideGrad
      iConBc(lbound)=1
    ENDIF

    CALL RFLO_SetIndexRange
    CALL RFLO_CopyPatchEdgeGrad

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

! ==============================================================================
! gradient completion subroutine for dummy cells at connecting boundaries
! ==============================================================================

CONTAINS

  SUBROUTINE RFLO_SetIndexRange

! set index-range be used in RFLO_CalcSideGrad and RFLO_CopyPatchEdgeGrad ----
      
    IF (lbound==1 .OR. lbound==3 .OR. lbound==5) THEN
      inbeg = ibeg-idir*(region%nDumCells-1)
      jnbeg = jbeg-jdir*(region%nDumCells-1)
      knbeg = kbeg-kdir*(region%nDumCells-1)
      inend = iend-2*idir+1
      jnend = jend-2*jdir+1
      knend = kend-2*kdir+1
    ELSE
      inbeg = ibeg-2*idir
      jnbeg = jbeg-2*jdir
      knbeg = kbeg-2*kdir
      inend = iend-idir*(region%nDumCells-1)+1
      jnend = jend-jdir*(region%nDumCells-1)+1
      knend = kend-kdir*(region%nDumCells-1)+1
    ENDIF

  END SUBROUTINE RFLO_SetIndexRange

  SUBROUTINE RFLO_CalcSideGrad

! ... local variables

  INTEGER :: inode,jnode,knode

! compute gradients at sides of current region having connecting bc --------

  inode = 0
  jnode = 0
  knode = 0
  IF( lbound==2 .OR. lbound==4 .OR. lbound==6 ) THEN
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

  IF (((knend>knbeg).AND.(jnend>jnbeg).AND.(inend>inbeg)).OR. &
      (lbound==1).OR.(lbound==3).OR.(lbound==5)) THEN

    DO k=knbeg+knode,knend
      DO j=jnbeg+jnode,jnend
        DO i=inbeg+inode,inend  

          ijkC    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
          im1jkC  = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
          ijm1kC  = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
          ijkm1C  = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
          im1jm1kC= IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
          im1jkm1C= IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
          ijm1km1C= IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
          ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          im1jkN  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
          ijm1kN  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
          ijkm1N  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

! - gradients at I-Face

! --- i-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,im1jkN))
          fny = .5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,im1jkN))
          fnz = .5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,im1jkN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
            tvar=var(l,im1jkC)
            gradi(lx,im1jkN) = gradi(lx,im1jkN) - fnx*tvar
            gradi(ly,im1jkN) = gradi(ly,im1jkN) - fny*tvar
            gradi(lz,im1jkN) = gradi(lz,im1jkN) - fnz*tvar

            gradi(lx,ijkN) = gradi(lx,ijkN) + fnx*tvar
            gradi(ly,ijkN) = gradi(ly,ijkN) + fny*tvar
            gradi(lz,ijkN) = gradi(lz,ijkN) + fnz*tvar
          END DO

! --- j-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,im1jkN))
          fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,im1jkN))
          fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,im1jkN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
!            fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                   var(l,im1jm1kC)+var(l,ijm1kC))

            fface(l) = ack(1,ijkN)*var(l,im1jm1kC)+ack(2,ijkN)*var(l,ijm1kC)+ &
                       ack(3,ijkN)*var(l,ijkC)+ack(4,ijkN)*var(l,im1jkC)
                                 
            gradi(lx,ijm1kN) = gradi(lx,ijm1kN) - fnx*fface(l)
            gradi(ly,ijm1kN) = gradi(ly,ijm1kN) - fny*fface(l)
            gradi(lz,ijm1kN) = gradi(lz,ijm1kN) - fnz*fface(l)

            gradi(lx,ijkN) = gradi(lx,ijkN) + fnx*fface(l)
            gradi(ly,ijkN) = gradi(ly,ijkN) + fny*fface(l)
            gradi(lz,ijkN) = gradi(lz,ijkN) + fnz*fface(l)
          END DO

! --- k-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(sk(XCOORD,im1jkN)+sk(XCOORD,ijkN))
          fny = .5_RFREAL*(sk(YCOORD,im1jkN)+sk(YCOORD,ijkN))
          fnz = .5_RFREAL*(sk(ZCOORD,im1jkN)+sk(ZCOORD,ijkN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
!            fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                   var(l,im1jkm1C)+var(l,ijkm1C))

            fface(l) = acj(1,ijkN)*var(l,im1jkm1C)+acj(2,ijkN)*var(l,im1jkC)+ &
                       acj(3,ijkN)*var(l,ijkC)+acj(4,ijkN)*var(l,ijkm1C)

            gradi(lx,ijkm1N) = gradi(lx,ijkm1N) - fnx*fface(l)
            gradi(ly,ijkm1N) = gradi(ly,ijkm1N) - fny*fface(l)
            gradi(lz,ijkm1N) = gradi(lz,ijkm1N) - fnz*fface(l)

            gradi(lx,ijkN) = gradi(lx,ijkN) + fnx*fface(l)
            gradi(ly,ijkN) = gradi(ly,ijkN) + fny*fface(l)
            gradi(lz,ijkN) = gradi(lz,ijkN) + fnz*fface(l)
          END DO

! - gradients at J-Face

! --- j-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijm1kN))
          fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijm1kN))
          fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijm1kN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
            tvar=var(l,ijm1kC)
            gradj(lx,ijm1kN) = gradj(lx,ijm1kN) - fnx*tvar
            gradj(ly,ijm1kN) = gradj(ly,ijm1kN) - fny*tvar
            gradj(lz,ijm1kN) = gradj(lz,ijm1kN) - fnz*tvar

            gradj(lx,ijkN) = gradj(lx,ijkN) + fnx*tvar
            gradj(ly,ijkN) = gradj(ly,ijkN) + fny*tvar
            gradj(lz,ijkN) = gradj(lz,ijkN) + fnz*tvar
          END DO

! --- i-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,ijm1kN))
          fny = .5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,ijm1kN))
          fnz = .5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,ijm1kN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
!            fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                   var(l,im1jm1kC)+var(l,ijm1kC))

            fface(l) = ack(1,ijkN)*var(l,im1jm1kC)+ack(2,ijkN)*var(l,ijm1kC)+ &
                       ack(3,ijkN)*var(l,ijkC)+ack(4,ijkN)*var(l,im1jkC)
                                 
            gradj(lx,im1jkN) = gradj(lx,im1jkN) - fnx*fface(l)
            gradj(ly,im1jkN) = gradj(ly,im1jkN) - fny*fface(l)
            gradj(lz,im1jkN) = gradj(lz,im1jkN) - fnz*fface(l)

            gradj(lx,ijkN) = gradj(lx,ijkN) + fnx*fface(l)
            gradj(ly,ijkN) = gradj(ly,ijkN) + fny*fface(l)
            gradj(lz,ijkN) = gradj(lz,ijkN) + fnz*fface(l)
          END DO

! --- k-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(sk(XCOORD,ijm1kN)+sk(XCOORD,ijkN))
          fny = .5_RFREAL*(sk(YCOORD,ijm1kN)+sk(YCOORD,ijkN))
          fnz = .5_RFREAL*(sk(ZCOORD,ijm1kN)+sk(ZCOORD,ijkN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
!            fface(l) = .25_RFREAL*(var(l,ijm1kC)+var(l,ijkC)+ &
!                                   var(l,ijm1km1C)+var(l,ijkm1C))

            fface(l) = aci(1,ijkN)*var(l,ijm1km1C)+aci(2,ijkN)*var(l,ijkm1C)+ &
                       aci(3,ijkN)*var(l,ijkC)+aci(4,ijkN)*var(l,ijm1kC)
                                 
            gradj(lx,ijkm1N) = gradj(lx,ijkm1N) - fnx*fface(l)
            gradj(ly,ijkm1N) = gradj(ly,ijkm1N) - fny*fface(l)
            gradj(lz,ijkm1N) = gradj(lz,ijkm1N) - fnz*fface(l)

            gradj(lx,ijkN) = gradj(lx,ijkN) + fnx*fface(l)
            gradj(ly,ijkN) = gradj(ly,ijkN) + fny*fface(l)
            gradj(lz,ijkN) = gradj(lz,ijkN) + fnz*fface(l)
          END DO

! - gradients at K-Face

! --- k-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(sk(XCOORD,ijkN)+sk(XCOORD,ijkm1N))
          fny = .5_RFREAL*(sk(YCOORD,ijkN)+sk(YCOORD,ijkm1N))
          fnz = .5_RFREAL*(sk(ZCOORD,ijkN)+sk(ZCOORD,ijkm1N))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
            tvar=var(l,ijkm1C)
            gradk(lx,ijkm1N) = gradk(lx,ijkm1N) - fnx*tvar
            gradk(ly,ijkm1N) = gradk(ly,ijkm1N) - fny*tvar
            gradk(lz,ijkm1N) = gradk(lz,ijkm1N) - fnz*tvar

            gradk(lx,ijkN) = gradk(lx,ijkN) + fnx*tvar
            gradk(ly,ijkN) = gradk(ly,ijkN) + fny*tvar
            gradk(lz,ijkN) = gradk(lz,ijkN) + fnz*tvar
          END DO

! --- j-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijkm1N))
          fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijkm1N))
          fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijkm1N))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
!            fface(l) = .25_RFREAL*(var(l,ijkm1C)+var(l,ijkC)+ &
!                                   var(l,ijm1km1C)+var(l,ijm1kC))

            fface(l) = aci(1,ijkN)*var(l,ijm1km1C)+aci(2,ijkN)*var(l,ijkm1C)+ &
                       aci(3,ijkN)*var(l,ijkC)+aci(4,ijkN)*var(l,ijm1kC)

            gradk(lx,ijm1kN) = gradk(lx,ijm1kN) - fnx*fface(l)
            gradk(ly,ijm1kN) = gradk(ly,ijm1kN) - fny*fface(l)
            gradk(lz,ijm1kN) = gradk(lz,ijm1kN) - fnz*fface(l)

            gradk(lx,ijkN) = gradk(lx,ijkN) + fnx*fface(l)
            gradk(ly,ijkN) = gradk(ly,ijkN) + fny*fface(l)
            gradk(lz,ijkN) = gradk(lz,ijkN) + fnz*fface(l)
          END DO

! --- i-face of auxiliary control volume (face vectors pointed inwards/negative)

          fnx = .5_RFREAL*(si(XCOORD,ijkm1N)+si(XCOORD,ijkN))
          fny = .5_RFREAL*(si(YCOORD,ijkm1N)+si(YCOORD,ijkN))
          fnz = .5_RFREAL*(si(ZCOORD,ijkm1N)+si(ZCOORD,ijkN))

          DO l=iBegV,iEndV
            lx=l+iBegG-iBegV 
            ly=lx+nvar
            lz=ly+nvar
!            fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                   var(l,im1jkm1C)+var(l,ijkm1C))

            fface(l) = acj(1,ijkN)*var(l,im1jkm1C)+acj(2,ijkN)*var(l,im1jkC)+ &
                       acj(3,ijkN)*var(l,ijkC)+acj(4,ijkN)*var(l,ijkm1C)

            gradk(lx,im1jkN) = gradk(lx,im1jkN) - fnx*fface(l)
            gradk(ly,im1jkN) = gradk(ly,im1jkN) - fny*fface(l)
            gradk(lz,im1jkN) = gradk(lz,im1jkN) - fnz*fface(l)

            gradk(lx,ijkN) = gradk(lx,ijkN) + fnx*fface(l)
            gradk(ly,ijkN) = gradk(ly,ijkN) + fny*fface(l)
            gradk(lz,ijkN) = gradk(lz,ijkN) + fnz*fface(l)
          END DO

        END DO ! i
      END DO   ! j
    END DO     ! k

  END IF

  IF (lbound==2) THEN

! - finish the last i-face (face vectors pointed inwards/negative)

    DO k=knbeg,knend-1 
      DO j=jnbeg,jnend-1  

        im1jkC  = IndIJK(inend  ,j  ,k  ,iCOff,ijCOff)
        ijkN    = IndIJK(inend+1,j  ,k  ,iNOff,ijNOff)
        im1jkN  = IndIJK(inend  ,j  ,k  ,iNOff,ijNOff)

        fnx = .5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,im1jkN))
        fny = .5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,im1jkN))
        fnz = .5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,im1jkN))
        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
          tvar=var(l,im1jkC)
          gradi(lx,im1jkN) = gradi(lx,im1jkN) - fnx*tvar
          gradi(ly,im1jkN) = gradi(ly,im1jkN) - fny*tvar
          gradi(lz,im1jkN) = gradi(lz,im1jkN) - fnz*tvar
        END DO

      END DO ! j
    END DO   ! k

  ELSEIF (lbound==4) THEN

! - finish the last j-face (face vectors pointed inwards/negative)

    DO k=knbeg,knend-1
      DO i=inbeg,inend-1  
 
        ijm1kC  = IndIJK(i  ,jnend  ,k  ,iCOff,ijCOff)
        ijkN    = IndIJK(i  ,jnend+1,k  ,iNOff,ijNOff)
        ijm1kN  = IndIJK(i  ,jnend  ,k  ,iNOff,ijNOff)

        fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijm1kN))
        fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijm1kN))
        fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijm1kN))
        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
          tvar=var(l,ijm1kC)
          gradj(lx,ijm1kN) = gradj(lx,ijm1kN) - fnx*tvar
          gradj(ly,ijm1kN) = gradj(ly,ijm1kN) - fny*tvar
          gradj(lz,ijm1kN) = gradj(lz,ijm1kN) - fnz*tvar
        END DO

      END DO ! i
    END DO   ! k

  ELSEIF (lbound==6) THEN

! - finish the last k-face (face vectors pointed inwards/negative)

    DO j=jnbeg,jnend-1 
      DO i=inbeg,inend-1  

        ijkm1C  = IndIJK(i  ,j  ,knend  ,iCOff,ijCOff)
        ijkN    = IndIJK(i  ,j  ,knend+1,iNOff,ijNOff)
        ijkm1N  = IndIJK(i  ,j  ,knend  ,iNOff,ijNOff)

        fnx = .5_RFREAL*(sk(XCOORD,ijkN)+sk(XCOORD,ijkm1N))
        fny = .5_RFREAL*(sk(YCOORD,ijkN)+sk(YCOORD,ijkm1N))
        fnz = .5_RFREAL*(sk(ZCOORD,ijkN)+sk(ZCOORD,ijkm1N))
        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
          tvar=var(l,ijkm1C)
          gradk(lx,ijkm1N) = gradk(lx,ijkm1N) - fnx*tvar
          gradk(ly,ijkm1N) = gradk(ly,ijkm1N) - fny*tvar
          gradk(lz,ijkm1N) = gradk(lz,ijkm1N) - fnz*tvar
        END DO

      END DO ! i
    END DO   ! j

  END IF

! finally, division by face averaged volume

  DO k=knbeg,knend
    DO j=jnbeg,jnend
      DO i=inbeg,inend  

        ijkC    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        im1jkC  = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijm1kC  = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkm1C  = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
       
        rvol = 2.0_RFREAL/(vol(im1jkC)+vol(ijkC))
        DO l=iBegG,iEndG
          gradi(l,ijkN) = gradi(l,ijkN)*rvol
        END DO

        rvol = 2.0_RFREAL/(vol(ijm1kC)+vol(ijkC))
        DO l=iBegG,iEndG
          gradj(l,ijkN) = gradj(l,ijkN)*rvol
        END DO

        rvol = 2.0_RFREAL/(vol(ijkm1C)+vol(ijkC))
        DO l=iBegG,iEndG
          gradk(l,ijkN) = gradk(l,ijkN)*rvol
        END DO

      END DO ! i
    END DO   ! j
  END DO     ! k

  END SUBROUTINE RFLO_CalcSideGrad

 
  SUBROUTINE RFLO_CopyPatchEdgeGrad

! ... local variables

  INTEGER :: ijkN1

! copy gradients to patch edges ---------------------------------------------

  IF ((lbound==1).OR.(lbound==2)) THEN 
    DO i=inbeg,inend
      DO k=knbeg,knend-1
        ijkN    = IndIJK(i   ,jnbeg  ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(ibeg,jnbeg  ,k,iNOff,ijNOff)
        gradj(:,ijkN)=gradj(:,ijkN1)
        ijkN    = IndIJK(i   ,jnend  ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(ibeg,jnend  ,k,iNOff,ijNOff)
        gradj(:,ijkN)=gradj(:,ijkN1)
      END DO  
      DO j=jnbeg,jnend-1
        ijkN    = IndIJK(i   ,j,knbeg  ,iNOff,ijNOff)
        ijkN1   = IndIJK(ibeg,j,knbeg  ,iNOff,ijNOff)
        gradk(:,ijkN)=gradk(:,ijkN1)
        ijkN    = IndIJK(i   ,j,knend  ,iNOff,ijNOff)
        ijkN1   = IndIJK(ibeg,j,knend  ,iNOff,ijNOff)
        gradk(:,ijkN)=gradk(:,ijkN1)
      END DO  
    END DO   ! i
  ELSEIF ((lbound==3).OR.(lbound==4)) THEN
    DO j=jnbeg,jnend
      DO k=knbeg,knend-1
        ijkN    = IndIJK(inbeg  ,j   ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(inbeg  ,jbeg,k,iNOff,ijNOff)
        gradi(:,ijkN)=gradi(:,ijkN1)
        ijkN    = IndIJK(inend  ,j   ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(inend  ,jbeg,k,iNOff,ijNOff)
        gradi(:,ijkN)=gradi(:,ijkN1)
      END DO  
      DO i=inbeg,inend-1
        ijkN    = IndIJK(i,j   ,knbeg  ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jbeg,knbeg  ,iNOff,ijNOff)
        gradk(:,ijkN)=gradk(:,ijkN1)
        ijkN    = IndIJK(i,j   ,knend  ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jbeg,knend  ,iNOff,ijNOff)
        gradk(:,ijkN)=gradk(:,ijkN1)
      END DO  
    END DO   ! j
  ELSEIF ((lbound==5).OR.(lbound==6)) THEN
    DO k=knbeg,knend
      DO i=inbeg,inend-1
        ijkN    = IndIJK(i,jnbeg  ,k   ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jnbeg  ,kbeg,iNOff,ijNOff)
        gradj(:,ijkN)=gradj(:,ijkN1)
        ijkN    = IndIJK(i,jnend  ,k   ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jnend  ,kbeg,iNOff,ijNOff)
        gradj(:,ijkN)=gradj(:,ijkN1)
      END DO  
      DO j=jnbeg,jnend-1
        ijkN    = IndIJK(inbeg  ,j,k   ,iNOff,ijNOff)
        ijkN1   = IndIJK(inbeg  ,j,kbeg,iNOff,ijNOff)
        gradi(:,ijkN)=gradi(:,ijkN1)
        ijkN    = IndIJK(inend  ,j,k   ,iNOff,ijNOff)
        ijkN1   = IndIJK(inend  ,j,kbeg,iNOff,ijNOff)
        gradi(:,ijkN)=gradi(:,ijkN1)
      END DO  
    END DO   ! k
  END IF     ! lbound 

  END SUBROUTINE RFLO_CopyPatchEdgeGrad

END SUBROUTINE RFLO_CalcGradConnBc

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_CalcGradConnBc.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:15  wasistho
! lower to upper case
!
! Revision 1.11  2004/08/03 22:48:43  wasistho
! changed cell2edge averaging to grid dependent avg
!
! Revision 1.10  2004/01/22 03:56:27  wasistho
! extrapolated to patch edges from different direction of side faces
!
! Revision 1.9  2003/12/07 17:48:39  jiao
! Chaged in order to compile with Intel compilers on NCSA machines
!
! Revision 1.8  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.4  2003/01/10 17:58:42  jblazek
! Added missing explicit interfaces.
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
! Revision 1.2  2002/07/19 23:43:00  wasistho
! made compliant with CODING RULE
!
!******************************************************************************








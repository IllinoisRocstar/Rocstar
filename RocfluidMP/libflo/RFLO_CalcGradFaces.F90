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
! Purpose: compute gradients of any vector or scalar at i, j and k faces
!          in one big loop
!
! Description: gradients are computed by second order finite volume method
!              using auxiliary control volume whose center is cell-face-center
!              df/dx = 1/V int(f.(nx,0,0))dS
!              df/dy = 1/V int(f.(0,ny,0))dS
!              df/dz = 1/V int(f.(0,0,nz))dS
!
! Input: region        = data of current region.
!        ilev          = data of current level.
!        iBegV, iEndV  = begin and end var index
!        iBegG, iEndG  = begin and end gradient index
!        var           = variables, the gradient of which to be determined.
!
! Output: gradi,j,k = gradients of 'var' at i,j,k-face.
!
! Notes: the computation covers region-interior up to -boundaries/sides
!        this method also applied at dummy faces of connecting sides
!        (subroutine RFLO_CalcGradConnBc)
!
!******************************************************************************
!
! $Id: RFLO_CalcGradFaces.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradFaces( region,ilev,iBegV,iEndV,iBegG,iEndG,var, &
                               gradi,gradj,gradk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: ilev, iBegV, iEndV, iBegG, iEndG
  REAL(RFREAL), POINTER :: var(:,:), gradi(:,:),gradj(:,:),gradk(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, lx, ly, lz

! ... local variables
  INTEGER :: inbeg, inend, jnbeg, jnend, knbeg, knend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: ijkC, im1jkC, ijm1kC, ijkm1C, im1jm1kC, im1jkm1C, ijm1km1C
  INTEGER :: ijkN, im1jkN, ijm1kN, ijkm1N
  INTEGER :: nvar

  REAL(RFREAL) :: rvol, fnx, fny, fnz, fface(iBegV:iEndV)
  REAL(RFREAL), POINTER :: aci(:,:), acj(:,:), ack(:,:)
  REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:), vol(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradFaces',&
  'RFLO_CalcGradFaces.F90' )

! get dimensions and pointers -------------------------------------------------

  nvar = iEndV - iBegV + 1

  aci => region%levels(ilev)%grid%c2eCoI
  acj => region%levels(ilev)%grid%c2eCoJ
  ack => region%levels(ilev)%grid%c2eCoK
  si  => region%levels(ilev)%grid%si
  sj  => region%levels(ilev)%grid%sj
  sk  => region%levels(ilev)%grid%sk
  vol => region%levels(ilev)%grid%vol

! get cell and node dimensions ------------------------------------------------

  CALL RFLO_GetDimensPhysNodes( region,ilev,inbeg,inend, &
                                jnbeg,jnend,knbeg,knend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

! initialize gradients --------------------------------------------------------

  gradi(iBegG:iEndG,:) = 0.0_RFREAL
  gradj(iBegG:iEndG,:) = 0.0_RFREAL
  gradk(iBegG:iEndG,:) = 0.0_RFREAL

! begin gradients computation -------------------------------------------------

  DO k=knbeg,knend+1 
    DO j=jnbeg,jnend+1 
      DO i=inbeg,inend+1  

        ijkC    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        im1jkC  = ijkC-1
        ijm1kC  = ijkC-iCOff
        ijkm1C  = ijkC-ijCOff
        im1jm1kC= im1jkC-iCOff
        im1jkm1C= im1jkC-ijCOff
        ijm1km1C= ijm1kC-ijCOff
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        im1jkN  = ijkN-1
        ijm1kN  = ijkN-iNOff
        ijkm1N  = ijkN-ijNOff
        
! - gradients at I-Face

! --- i-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,im1jkN))
        fny = .5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,im1jkN))
        fnz = .5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,im1jkN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
          gradi(lx,im1jkN) = gradi(lx,im1jkN) - fnx*var(l,im1jkC)
          gradi(ly,im1jkN) = gradi(ly,im1jkN) - fny*var(l,im1jkC)
          gradi(lz,im1jkN) = gradi(lz,im1jkN) - fnz*var(l,im1jkC)

          gradi(lx,ijkN) = gradi(lx,ijkN) + fnx*var(l,im1jkC)
          gradi(ly,ijkN) = gradi(ly,ijkN) + fny*var(l,im1jkC)
          gradi(lz,ijkN) = gradi(lz,ijkN) + fnz*var(l,im1jkC)
        ENDDO

! --- j-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,im1jkN))
        fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,im1jkN))
        fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,im1jkN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
!          fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                 var(l,im1jm1kC)+var(l,ijm1kC))

          fface(l) = ack(1,ijkN)*var(l,im1jm1kC)+ack(2,ijkN)*var(l,ijm1kC)+ &
                     ack(3,ijkN)*var(l,ijkC)+ack(4,ijkN)*var(l,im1jkC)
                                 
          gradi(lx,ijm1kN) = gradi(lx,ijm1kN) - fnx*fface(l)
          gradi(ly,ijm1kN) = gradi(ly,ijm1kN) - fny*fface(l)
          gradi(lz,ijm1kN) = gradi(lz,ijm1kN) - fnz*fface(l)

          gradi(lx,ijkN) = gradi(lx,ijkN) + fnx*fface(l)
          gradi(ly,ijkN) = gradi(ly,ijkN) + fny*fface(l)
          gradi(lz,ijkN) = gradi(lz,ijkN) + fnz*fface(l)
        ENDDO

! --- k-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(sk(XCOORD,im1jkN)+sk(XCOORD,ijkN))
        fny = .5_RFREAL*(sk(YCOORD,im1jkN)+sk(YCOORD,ijkN))
        fnz = .5_RFREAL*(sk(ZCOORD,im1jkN)+sk(ZCOORD,ijkN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
!          fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                 var(l,im1jkm1C)+var(l,ijkm1C))

          fface(l) = acj(1,ijkN)*var(l,im1jkm1C)+acj(2,ijkN)*var(l,im1jkC)+ &
                     acj(3,ijkN)*var(l,ijkC)+acj(4,ijkN)*var(l,ijkm1C)

          gradi(lx,ijkm1N) = gradi(lx,ijkm1N) - fnx*fface(l)
          gradi(ly,ijkm1N) = gradi(ly,ijkm1N) - fny*fface(l)
          gradi(lz,ijkm1N) = gradi(lz,ijkm1N) - fnz*fface(l)

          gradi(lx,ijkN) = gradi(lx,ijkN) + fnx*fface(l)
          gradi(ly,ijkN) = gradi(ly,ijkN) + fny*fface(l)
          gradi(lz,ijkN) = gradi(lz,ijkN) + fnz*fface(l)
        ENDDO
        
! - gradients at J-Face

! --- j-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijm1kN))
        fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijm1kN))
        fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijm1kN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
          gradj(lx,ijm1kN) = gradj(lx,ijm1kN) - fnx*var(l,ijm1kC)
          gradj(ly,ijm1kN) = gradj(ly,ijm1kN) - fny*var(l,ijm1kC)
          gradj(lz,ijm1kN) = gradj(lz,ijm1kN) - fnz*var(l,ijm1kC)

          gradj(lx,ijkN) = gradj(lx,ijkN) + fnx*var(l,ijm1kC)
          gradj(ly,ijkN) = gradj(ly,ijkN) + fny*var(l,ijm1kC)
          gradj(lz,ijkN) = gradj(lz,ijkN) + fnz*var(l,ijm1kC)
        ENDDO

! --- i-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,ijm1kN))
        fny = .5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,ijm1kN))
        fnz = .5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,ijm1kN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
!          fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                 var(l,im1jm1kC)+var(l,ijm1kC))

          fface(l) = ack(1,ijkN)*var(l,im1jm1kC)+ack(2,ijkN)*var(l,ijm1kC)+ &
                     ack(3,ijkN)*var(l,ijkC)+ack(4,ijkN)*var(l,im1jkC)
                                 
          gradj(lx,im1jkN) = gradj(lx,im1jkN) - fnx*fface(l)
          gradj(ly,im1jkN) = gradj(ly,im1jkN) - fny*fface(l)
          gradj(lz,im1jkN) = gradj(lz,im1jkN) - fnz*fface(l)

          gradj(lx,ijkN) = gradj(lx,ijkN) + fnx*fface(l)
          gradj(ly,ijkN) = gradj(ly,ijkN) + fny*fface(l)
          gradj(lz,ijkN) = gradj(lz,ijkN) + fnz*fface(l)
        ENDDO

! --- k-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(sk(XCOORD,ijm1kN)+sk(XCOORD,ijkN))
        fny = .5_RFREAL*(sk(YCOORD,ijm1kN)+sk(YCOORD,ijkN))
        fnz = .5_RFREAL*(sk(ZCOORD,ijm1kN)+sk(ZCOORD,ijkN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
!          fface(l) = .25_RFREAL*(var(l,ijm1kC)+var(l,ijkC)+ &
!                                 var(l,ijm1km1C)+var(l,ijkm1C))

          fface(l) = aci(1,ijkN)*var(l,ijm1km1C)+aci(2,ijkN)*var(l,ijkm1C)+ &
                     aci(3,ijkN)*var(l,ijkC)+aci(4,ijkN)*var(l,ijm1kC)
                                 
          gradj(lx,ijkm1N) = gradj(lx,ijkm1N) - fnx*fface(l)
          gradj(ly,ijkm1N) = gradj(ly,ijkm1N) - fny*fface(l)
          gradj(lz,ijkm1N) = gradj(lz,ijkm1N) - fnz*fface(l)

          gradj(lx,ijkN) = gradj(lx,ijkN) + fnx*fface(l)
          gradj(ly,ijkN) = gradj(ly,ijkN) + fny*fface(l)
          gradj(lz,ijkN) = gradj(lz,ijkN) + fnz*fface(l)
        ENDDO
        
! - gradients at K-Face

! --- k-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(sk(XCOORD,ijkN)+sk(XCOORD,ijkm1N))
        fny = .5_RFREAL*(sk(YCOORD,ijkN)+sk(YCOORD,ijkm1N))
        fnz = .5_RFREAL*(sk(ZCOORD,ijkN)+sk(ZCOORD,ijkm1N))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
          gradk(lx,ijkm1N) = gradk(lx,ijkm1N) - fnx*var(l,ijkm1C)
          gradk(ly,ijkm1N) = gradk(ly,ijkm1N) - fny*var(l,ijkm1C)
          gradk(lz,ijkm1N) = gradk(lz,ijkm1N) - fnz*var(l,ijkm1C)

          gradk(lx,ijkN) = gradk(lx,ijkN) + fnx*var(l,ijkm1C)
          gradk(ly,ijkN) = gradk(ly,ijkN) + fny*var(l,ijkm1C)
          gradk(lz,ijkN) = gradk(lz,ijkN) + fnz*var(l,ijkm1C)
        ENDDO

! --- j-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijkm1N))
        fny = .5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijkm1N))
        fnz = .5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijkm1N))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
!          fface(l) = .25_RFREAL*(var(l,ijkm1C)+var(l,ijkC)+ &
!                                 var(l,ijm1km1C)+var(l,ijm1kC))

          fface(l) = aci(1,ijkN)*var(l,ijm1km1C)+aci(2,ijkN)*var(l,ijkm1C)+ &
                     aci(3,ijkN)*var(l,ijkC)+aci(4,ijkN)*var(l,ijm1kC)

          gradk(lx,ijm1kN) = gradk(lx,ijm1kN) - fnx*fface(l)
          gradk(ly,ijm1kN) = gradk(ly,ijm1kN) - fny*fface(l)
          gradk(lz,ijm1kN) = gradk(lz,ijm1kN) - fnz*fface(l)

          gradk(lx,ijkN) = gradk(lx,ijkN) + fnx*fface(l)
          gradk(ly,ijkN) = gradk(ly,ijkN) + fny*fface(l)
          gradk(lz,ijkN) = gradk(lz,ijkN) + fnz*fface(l)
        ENDDO

! --- i-face of auxiliary control volume (face vectors pointed inwards/negative)

        fnx = .5_RFREAL*(si(XCOORD,ijkm1N)+si(XCOORD,ijkN))
        fny = .5_RFREAL*(si(YCOORD,ijkm1N)+si(YCOORD,ijkN))
        fnz = .5_RFREAL*(si(ZCOORD,ijkm1N)+si(ZCOORD,ijkN))

        DO l=iBegV,iEndV
          lx=l+iBegG-iBegV 
          ly=lx+nvar
          lz=ly+nvar
!          fface(l) = .25_RFREAL*(var(l,im1jkC)+var(l,ijkC)+ &
!                                 var(l,im1jkm1C)+var(l,ijkm1C))

          fface(l) = acj(1,ijkN)*var(l,im1jkm1C)+acj(2,ijkN)*var(l,im1jkC)+ &
                     acj(3,ijkN)*var(l,ijkC)+acj(4,ijkN)*var(l,ijkm1C)

          gradk(lx,im1jkN) = gradk(lx,im1jkN) - fnx*fface(l)
          gradk(ly,im1jkN) = gradk(ly,im1jkN) - fny*fface(l)
          gradk(lz,im1jkN) = gradk(lz,im1jkN) - fnz*fface(l)

          gradk(lx,ijkN) = gradk(lx,ijkN) + fnx*fface(l)
          gradk(ly,ijkN) = gradk(ly,ijkN) + fny*fface(l)
          gradk(lz,ijkN) = gradk(lz,ijkN) + fnz*fface(l)
        ENDDO

      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! - division by face volume

  DO k=knbeg,knend
    DO j=jnbeg,jnend
      DO i=inbeg,inend  

        ijkC    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        im1jkC  = ijkC-1
        ijm1kC  = ijkC-iCOff
        ijkm1C  = ijkC-ijCOff
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
       
        rvol = 2.0_RFREAL/(vol(im1jkC)+vol(ijkC))
        DO l=iBegG,iEndG
          gradi(l,ijkN) = gradi(l,ijkN)*rvol
        ENDDO

        rvol = 2.0_RFREAL/(vol(ijm1kC)+vol(ijkC))
        DO l=iBegG,iEndG
          gradj(l,ijkN) = gradj(l,ijkN)*rvol
        ENDDO

        rvol = 2.0_RFREAL/(vol(ijkm1C)+vol(ijkC))
        DO l=iBegG,iEndG
          gradk(l,ijkN) = gradk(l,ijkN)*rvol
        ENDDO

      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! provide approximate gradients at the edge boundaries by copying
! from the adjacent cells

  CALL RFLO_CopyEdgeGrad

! up to this point gradi,j,k are well defined in the following range for 
! block interfaces and periodic bc: [knbeg:knend,jnbeg:jnend,inbeg:inend];
! there are incomplete contributions from the auxiliary control volumes
! at the boundaries which will be completed in the dummy treatment of the
! connecting boundary conditions;
! for physical bc, special treatment are needed for the gradients at the patches,
! depending on whether dummy variables are linear or  zeroth extrapolated 
! (in slip_wall and injection bc).

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

!===================================================================
! subroutine providing approximate gradients at the edge boundaries 
!===================================================================

CONTAINS
 
  SUBROUTINE RFLO_CopyEdgeGrad

! ... local variables

  INTEGER :: ijkN1

! gradients at juctions edges and corners approximated by copying from
! adjacent interior faces

! lbound=2
  DO i=inend,inend
    DO k=knbeg,knend-1
       ijkN    = IndIJK(i  ,jnbeg  ,k,iNOff,ijNOff)
       ijkN1   = ijkN+iNOff
       gradj(:,ijkN)=gradj(:,ijkN1)
       ijkN    = IndIJK(i  ,jnend  ,k,iNOff,ijNOff)
       ijkN1   = ijkN-iNOff
       gradj(:,ijkN)=gradj(:,ijkN1)
    END DO  
    DO j=jnbeg,jnend-1
       ijkN    = IndIJK(i  ,j,knbeg  ,iNOff,ijNOff)
       ijkN1   = ijkN+ijNOff
       gradk(:,ijkN)=gradk(:,ijkN1)
       ijkN    = IndIJK(i  ,j,knend  ,iNOff,ijNOff)
       ijkN1   = ijkN-ijNOff
       gradk(:,ijkN)=gradk(:,ijkN1)
    ENDDO  
  ENDDO

! lbound=4
  DO j=jnend,jnend
    DO k=knbeg,knend-1
       ijkN    = IndIJK(inbeg  ,j  ,k,iNOff,ijNOff)
       ijkN1   = ijkN+1
       gradi(:,ijkN)=gradi(:,ijkN1)
       ijkN    = IndIJK(inend  ,j  ,k,iNOff,ijNOff)
       ijkN1   = ijkN-1
       gradi(:,ijkN)=gradi(:,ijkN1)
    ENDDO  
    DO i=inbeg,inend-1
       ijkN    = IndIJK(i,j  ,knbeg  ,iNOff,ijNOff)
       ijkN1   = ijkN+ijNOff
       gradk(:,ijkN)=gradk(:,ijkN1)
       ijkN    = IndIJK(i,j  ,knend  ,iNOff,ijNOff)
       ijkN1   = ijkN-ijNOff
       gradk(:,ijkN)=gradk(:,ijkN1)
    ENDDO  
  ENDDO  

! lbound=6
  DO k=knend,knend
    DO i=inbeg,inend-1
       ijkN    = IndIJK(i,jnbeg  ,k  ,iNOff,ijNOff)
       ijkN1   = ijkN+iNOff
       gradj(:,ijkN)=gradj(:,ijkN1)
       ijkN    = IndIJK(i,jnend  ,k  ,iNOff,ijNOff)
       ijkN1   = ijkN-iNOff
       gradj(:,ijkN)=gradj(:,ijkN1)
    ENDDO  
    DO j=jnbeg,jnend-1
       ijkN    = IndIJK(inbeg  ,j,k  ,iNOff,ijNOff)
       ijkN1   = ijkN+1
       gradi(:,ijkN)=gradi(:,ijkN1)
       ijkN    = IndIJK(inend  ,j,k  ,iNOff,ijNOff)
       ijkN1   = ijkN-1
       gradi(:,ijkN)=gradi(:,ijkN1)
    ENDDO  
  ENDDO  

  END SUBROUTINE RFLO_CopyEdgeGrad

END SUBROUTINE RFLO_CalcGradFaces

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_CalcGradFaces.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:13  haselbac
! Removed tabs
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.10  2004/08/03 22:48:36  wasistho
! changed cell2edge averaging to grid dependent avg
!
! Revision 1.9  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/07/31 02:05:50  wasistho
!  put bound iBegG:iEndG in zero initialisation
!
! Revision 1.5  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.4  2003/01/10 17:58:42  jblazek
! Added missing explicit interfaces.
!
! Revision 1.3  2002/10/14 23:47:25  wasistho
! Minor tunning to get speedup
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/09/02 22:58:54  wasistho
! RFLO grad routines migrated from rocflo to libflo
!
! Revision 1.4  2002/08/29 18:36:15  wasistho
! Generalized range of fface
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:44:12  wasistho
! made compliant with CODING RULE
!
!******************************************************************************








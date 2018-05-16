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
! Purpose: Compute gradients of a given vector at i,j,k faces
!
! Description: Computation performed first at interior points 
!              then at patches for:
!              - connecting (periodic and block interface) b.c. case 
!              - physical b.c. case 
!              then their corresponding dummies
!
! Input: region = data of current region
!        iBegV  = begin component index of vector variable to be grad
!        iEndV  = end   component index of vector variable to be grad
!        iBegG  = begin component index of resulting gradient vector
!        iEndG  = end   component index of resulting gradient vector
!        var    = vector variable to be grad
!
! Output: gradi,gradj,gradk = region%levels%mixt%gradi,j,k =
!         resulting gradient vector at i,j,k faces
!
! Notes: Subroutine CheckVelTIndx and CheckGrad contained in this file 
!        to check consistency of component indices and to check results of 
!        velocity and temperature gradients-computation
!
!******************************************************************************
!
! $Id: RFLO_CalcGradVector.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradVector( region,iBegV,iEndV,iBegG,iEndG,var, &
                                gradi,gradj,gradk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_CalcGradFaces, &
        RFLO_CalcGradConnBc, RFLO_CalcGradPhysBc, RFLO_CalcGradDummy, &
        RFLO_CopyEdgeFaceNorm, RFLO_CopyEdgeFaceParal 
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: iBegV, iEndV, iBegG, iEndG
  REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:),  gradk(:,:)

! ... loop variables
  INTEGER :: i, j, k, ipatch

! ... local variables
  INTEGER :: ilev
  INTEGER :: iConBc(6)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradVector',&
  'RFLO_CalcGradVector.F90' )

! get some parameters and flags --------------------------------------------

  ilev = region%currLevel
  iConBc(:) = 0   ! this flag identifies connecting patch is treated or not yet

! check consistency of indices

  CALL CheckIndxRange

! gradient in i, j and k-faces without dummy faces -------------------------

  CALL RFLO_CalcGradFaces( region,ilev,iBegV,iEndV,iBegG,iEndG,var, &
                           gradi,gradj,gradk )

! apply boundary conditions to the gradients at the patches -----------------
! ConnBc covers connecting bc: periodic and region-interface
! PhysBc covers physical bc: in/outflow, farfield, wall and injection
! Symmetry bc is covered in PhysBc for compactness 

  DO ipatch=1,region%nPatches
    CALL RFLO_CalcGradConnBc( region,region%levels(iLev)%patches(ipatch), &
                              iConBc,iBegV,iEndV,iBegG,iEndG,var, &
                              gradi,gradj,gradk )
    CALL RFLO_CalcGradPhysBc( region,region%levels(iLev)%patches(ipatch), & 
                              iBegV,iEndV,iBegG,iEndG,var,gradi,gradj,gradk )
  ENDDO

! extend to patch dummies

  DO ipatch=1,region%nPatches
    CALL RFLO_CalcGradDummy( region,region%levels(iLev)%patches(ipatch), & 
                             iBegV,iEndV,iBegG,iEndG,gradi,gradj,gradk )
  ENDDO

! and copy to edge/corner dummies -------------------------------------------- 
! the edge length is extended that corners are covered to minimize work

  CALL RFLO_CopyEdgeFaceNorm( region,iBegG,iEndG,gradi,gradj,gradk )
  CALL RFLO_CopyEdgeFaceParal( region,iBegG,iEndG,gradi,gradj,gradk )

! if needed, check all gradients up to dummy points (for debugging)-----------

#ifdef CHECK_GRAD
  CALL CheckGrad
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

!==============================================================================

CONTAINS

  SUBROUTINE CheckIndxRange
    USE ModError

! ... local variables
    INTEGER :: nvar, ngrd, istop

! - get number of array components and check consistency of component index

    nvar = iEndV - iBegV+1
    ngrd = iEndG - iBegG+1

    istop = 0
    IF ((nvar<1) .OR. (ngrd<1) .OR. (MOD(ngrd,nvar)/=0)) istop = 1
    IF (istop == 1) THEN 
      CALL ErrorStop( region%global,ERR_GRAD_INDEX,__LINE__, &
                      'Gradient indices are not consistent' )
    ENDIF

  END SUBROUTINE CheckIndxRange

!==============================================================================

#ifdef CHECK_GRAD
  SUBROUTINE CheckGrad

    USE ModGlobal, ONLY     : t_global
    USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
                              RFLO_GetCellOffset
    USE ModParameters
    IMPLICIT NONE

#include "Indexing.h"

    INTEGER      :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
    INTEGER      :: ijkC,ijkN,iCOff,iNOff,ijCOff,ijNOff
    INTEGER      :: ip1jkN,ijp1kN,ijkp1N,ip1jp1kN,ip1jkp1N,ijp1kp1N,ip1jp1kp1N
    INTEGER      :: m,l, mindx(9),iermax(9,3),jermax(9,3),kermax(9,3)
    REAL(RFREAL) :: xc, yc, zc
    REAL(RFREAL) :: diff(9,3), errmax(9,3), rval(9)

    CALL RFLO_GetDimensDummyNodes( region,ilev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )
    CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

    mindx(1)=1  ;rval(1)=0._RFREAL
    mindx(2)=2  ;rval(2)=0._RFREAL
    mindx(3)=3  ;rval(3)=0._RFREAL 
    mindx(4)=5  ;rval(4)=0._RFREAL
    mindx(5)=6  ;rval(5)=0._RFREAL
    mindx(6)=7  ;rval(6)=0._RFREAL
    mindx(7)=9  ;rval(7)=0._RFREAL
    mindx(8)=10 ;rval(8)=0._RFREAL
    mindx(9)=11 ;rval(9)=0._RFREAL
    DO m=1,9
      errmax(m,1)=0.0_RFREAL
      errmax(m,2)=0.0_RFREAL
      errmax(m,3)=0.0_RFREAL
    ENDDO
    DO k=kdnbeg,kdnend; DO j=jdnbeg,jdnend; DO i=idnbeg,idnend;  
      ijkN      = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
      ip1jkN    = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
      ijp1kN    = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
      ijkp1N    = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
      ip1jp1kN  = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
      ip1jkp1N  = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
      ijp1kp1N  = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
      ip1jp1kp1N= IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
      ijkC      = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
      xc=(region%levels(iLev)%grid%xyz(XCOORD,ijkN)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ip1jkN)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ijp1kN)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ijkp1N)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ip1jp1kN)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ijp1kp1N)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ip1jkp1N)+ &
          region%levels(iLev)%grid%xyz(XCOORD,ip1jp1kp1N))/8._RFREAL
      yc=(region%levels(iLev)%grid%xyz(YCOORD,ijkN)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ip1jkN)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ijp1kN)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ijkp1N)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ip1jp1kN)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ijp1kp1N)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ip1jkp1N)+ &
          region%levels(iLev)%grid%xyz(YCOORD,ip1jp1kp1N))/8._RFREAL
      zc=(region%levels(iLev)%grid%xyz(ZCOORD,ijkN)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ip1jkN)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ijp1kN)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ijkp1N)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ip1jp1kN)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ijp1kp1N)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ip1jkp1N)+ &
          region%levels(iLev)%grid%xyz(ZCOORD,ip1jp1kp1N))/8._RFREAL

      DO m=1,9
! ----- averaged cell values
!        diff(m,1)=0.5d0*(gradi(mindx(m),ijkN)+gradi(mindx(m),ip1jkN))-rval(m)
!        diff(m,2)=0.5d0*(gradj(mindx(m),ijkN)+gradj(mindx(m),ijp1kN))-rval(m)
!        diff(m,3)=0.5d0*(gradk(mindx(m),ijkN)+gradk(mindx(m),ijkp1N))-rval(m)
! ----- face values
        diff(m,1) = gradi(mindx(m),ijkN)-rval(m)
        diff(m,2) = gradj(mindx(m),ijkN)-rval(m)
        diff(m,3) = gradk(mindx(m),ijkN)-rval(m)

        IF (ABS(diff(m,1)) > errmax(m,1)) THEN
          errmax(m,1) = ABS(diff(m,1))
          iermax(m,1) = i
          jermax(m,1) = j
          kermax(m,1) = k
        ENDIF
        IF (ABS(diff(m,2)) > errmax(m,2)) THEN
          errmax(m,2) = ABS(diff(m,2))
          iermax(m,2) = i
          jermax(m,2) = j
          kermax(m,2) = k
        ENDIF
        IF (ABS(diff(m,3)) > errmax(m,3)) THEN
          errmax(m,3) = ABS(diff(m,3))
          iermax(m,3) = i
          jermax(m,3) = j
          kermax(m,3) = k
        ENDIF
      ENDDO   ! m

      WRITE(STDOUT,*) i,j,k,gradi(1,ijkN),gradi(2,ijkN),gradi(3,ijkN)
      WRITE(STDOUT,*) i,j,k,gradi(5,ijkN),gradi(6,ijkN),gradi(7,ijkN)
      WRITE(STDOUT,*) i,j,k,gradi(9,ijkN),gradi(10,ijkN),gradi(11,ijkN)

      WRITE(STDOUT,*) i,j,k,gradj(1,ijkN),gradj(2,ijkN),gradj(3,ijkN)
      WRITE(STDOUT,*) i,j,k,gradj(5,ijkN),gradj(6,ijkN),gradj(7,ijkN)
      WRITE(STDOUT,*) i,j,k,gradj(9,ijkN),gradj(10,ijkN),gradj(11,ijkN)

      WRITE(STDOUT,*) i,j,k,gradk(1,ijkN),gradk(2,ijkN),gradk(3,ijkN)
      WRITE(STDOUT,*) i,j,k,gradk(5,ijkN),gradk(6,ijkN),gradk(7,ijkN)
      WRITE(STDOUT,*) i,j,k,gradk(9,ijkN),gradk(10,ijkN),gradk(11,ijkN)
    ENDDO; ENDDO; ENDDO;

    WRITE(STDOUT,*)
    WRITE(STDOUT,*) 'region:',region%localNumber
    DO l=1,3
      WRITE(STDOUT,*)
      WRITE(STDOUT,*) 'max error in grad',l,'(ux) in i,j,k = ', &
                      errmax(1,l),iermax(1,l),jermax(1,l),kermax(1,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(vx) in i,j,k = ', &
                      errmax(2,l),iermax(2,l),jermax(2,l),kermax(2,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(wx) in i,j,k = ', &
                      errmax(3,l),iermax(3,l),jermax(3,l),kermax(3,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(uy) in i,j,k = ', &
                      errmax(4,l),iermax(4,l),jermax(4,l),kermax(4,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(vy) in i,j,k = ', &
                      errmax(5,l),iermax(5,l),jermax(5,l),kermax(5,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(wy) in i,j,k = ', &
                      errmax(6,l),iermax(6,l),jermax(6,l),kermax(6,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(uz) in i,j,k = ', &
                      errmax(7,l),iermax(7,l),jermax(7,l),kermax(7,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(vz) in i,j,k = ', &
                      errmax(8,l),iermax(8,l),jermax(8,l),kermax(8,l)
      WRITE(STDOUT,*) 'max error in grad',l,'(wz) in i,j,k = ', &
                      errmax(9,l),iermax(9,l),jermax(9,l),kermax(9,l)
    ENDDO

    IF (region%localNumber == region%global%nRegions) STOP

  END SUBROUTINE CheckGrad
#endif

END SUBROUTINE RFLO_CalcGradVector

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_CalcGradVector.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.10  2004/05/03 15:07:45  jferry
! minor bug fix:  changed "global" to "t_global"
!
! Revision 1.9  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.5  2003/01/10 17:58:42  jblazek
! Added missing explicit interfaces.
!
! Revision 1.4  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/09/03 00:00:52  wasistho
! Renamed CheckVelTIndx to CheckIndxRange
!
! Revision 1.2  2002/09/02 23:33:02  wasistho
! Moved grad index check to rocflo/RFLO_checkUserInput
!
! Revision 1.1  2002/09/02 22:58:54  wasistho
! RFLO grad routines migrated from rocflo to libflo
!
! Revision 1.2  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.1  2002/08/01 01:46:35  wasistho
! Renamed from libfloflu/calcGradVelT
!
!******************************************************************************








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
! Purpose: Initialisation of channel flow solution
!
! Description: initial solution is laminar channel flow
!
! Input: region = data of current region
!
! Output: channel solution to start simulation
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_coCnlInitSolutionFlo.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CoCnlInitSolution( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : MixtPerf_Eo_DGPVm
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE (t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, n

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER  :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER  :: iCOff, ijCOff, iNOff, ijNOff, ijkC, ijkCr, ijkN, ijkN1
  INTEGER  :: iCorner(8)
  REAL(RFREAL), POINTER :: xyz(:,:), cv(:,:), dv(:,:)
  REAL(RFREAL) :: rgas, yp, ym, yc, dy, muel, refMeanPgrad, delta
  REAL(RFREAL) :: rho, pres, Vm, Eo, xyzCell(3), twopi
  REAL(RFREAL) :: lambda,ampl,umax,nXwav,nYwav,nZwav,span,extn,fy,yscale,ywaves

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coCnlInitSolutionFlo.F90,v $ $Revision: 1.4 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_CoCnlInitSolution',&
  'PERI_coCnlInitSolutionFlo.F90' )

! get parameters and pointers ----------------------------------------------

  iLev   =  region%currLevel
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz    => region%levels(iLev)%grid%xyz
  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv

  IF ((global%flowType==FLOW_UNSTEADY .AND. &
       global%currentTime < PERI_REAL_SMALL) .OR. &
      (global%flowType==FLOW_STEADY .AND. global%currentIter==0)) THEN

! - get parameters and set conservative variables
    
    delta = global%refLength
    muel  = global%refVisc
    rho   = global%refDensity
    pres  = global%refPressure
    refMeanPgrad = region%periInput%meanPgrad
    cv(CV_MIXT_DENS,:) = rho
    cv(CV_MIXT_YMOM,:) = 0._RFREAL
    cv(CV_MIXT_ZMOM,:) = 0._RFREAL

    DO k = kpcbeg, kpcbeg
      DO i = ipcbeg, ipcbeg
        DO j = jpcbeg-region%nDumCells, jpcend+region%nDumCells
          ijkC  = IndIJK( i , j  , k ,iCOff,ijCOff) 
          ijkN  = IndIJK( i , j  , k ,iNOff,ijNOff) 
          ijkN1 = IndIJK( i , j+1, k ,iNOff,ijNOff) 
          dy    = xyz(YCOORD,ijkN1)-xyz(YCOORD,ijkN)
          yp    = xyz(YCOORD,ijkN1)
          ym    = xyz(YCOORD,ijkN) 
          yc    = 0.5_RFREAL*(xyz(YCOORD,ijkN1)+xyz(YCOORD,ijkN))
          IF (dy < 0._RFREAL) THEN
            CALL ErrorStop( region%global,ERR_PERI_GEO,__LINE__, &
                            'grid contains negative spacing' )
          ENDIF

          dv(DV_MIXT_UVEL,ijkC) = -0.5_RFREAL/muel*refMeanPgrad* &
                                   (delta*delta - yc*yc)

          cv(CV_MIXT_XMOM,ijkC) = dv(DV_MIXT_UVEL,ijkC)*rho

          Vm = SQRT( cv(CV_MIXT_XMOM,ijkC)*cv(CV_MIXT_XMOM,ijkC) + &
                     cv(CV_MIXT_YMOM,ijkC)*cv(CV_MIXT_YMOM,ijkC) + &	
                     cv(CV_MIXT_ZMOM,ijkC)*cv(CV_MIXT_ZMOM,ijkC))/rho
          Eo = MixtPerf_Eo_DGPVm(rho,global%refGamma,pres,Vm)
          cv(CV_MIXT_ENER,ijkC) = rho*Eo
        END DO
        DO j = 1, region%nDumCells 
          ijkN  = IndIJK( i , jpcbeg    , k ,iNOff,ijNOff) 
          ijkN1 = IndIJK( i , jpcend+1  , k ,iNOff,ijNOff) 
          yp    = xyz(YCOORD,ijkN1)
          ym    = xyz(YCOORD,ijkN) 
          ijkCr = IndIJK( i , jpcbeg+j-1, k ,iCOff,ijCOff) 
          ijkC  = IndIJK( i , jpcbeg-j  , k ,iCOff,ijCOff) 
          IF (ym == region%periInput%minmax(1)) THEN
            cv(CV_MIXT_XMOM,ijkC) = -cv(CV_MIXT_XMOM,ijkCr)
          ENDIF
          ijkCr = IndIJK( i , jpcend-j+1, k ,iCOff,ijCOff) 
          ijkC  = IndIJK( i , jpcend+j  , k ,iCOff,ijCOff) 
          IF (yp == region%periInput%minmax(2)) THEN
            cv(CV_MIXT_XMOM,ijkC) = -cv(CV_MIXT_XMOM,ijkCr)
          ENDIF
        END DO
      END DO
    END DO 

    DO k = kpcbeg-region%nDumCells, kpcend+region%nDumCells
      DO j = jpcbeg-region%nDumCells, jpcend+region%nDumCells 
        DO i = ipcbeg-region%nDumCells, ipcend+region%nDumCells 
          ijkC  = IndIJK( i    , j , k     ,iCOff,ijCOff) 
          ijkCr = IndIJK(ipcbeg, j ,kpcbeg ,iCOff,ijCOff) 
          cv(CV_MIXT_DENS,ijkC) = cv(CV_MIXT_DENS,ijkCr)
          cv(CV_MIXT_XMOM,ijkC) = cv(CV_MIXT_XMOM,ijkCr)
          cv(CV_MIXT_YMOM,ijkC) = cv(CV_MIXT_YMOM,ijkCr)
          cv(CV_MIXT_ZMOM,ijkC) = cv(CV_MIXT_ZMOM,ijkCr)
          cv(CV_MIXT_ENER,ijkC) = cv(CV_MIXT_ENER,ijkCr)
        END DO
      END DO
    END DO 

! - add perturbations

    lambda = 5._RFREAL
    ampl   = 0.3_RFREAL
    umax   = -0.5_RFREAL/muel*refMeanPgrad*delta*delta
    nXwav  = 3._RFREAL
    nYwav  = 2._RFREAL
    nZwav  = 2._RFREAL
    span   = xyz(ZCOORD,IndIJK( ipcbeg,jpcbeg,kpcend+1,iNOff,ijNOff ))- & 
             xyz(ZCOORD,IndIJK( ipcbeg,jpcbeg,kpcbeg  ,iNOff,ijNOff ))
    extn   = xyz(XCOORD,IndIJK( ipcend+1,jpcbeg,kpcbeg,iNOff,ijNOff ))- & 
             xyz(XCOORD,IndIJK( ipcbeg  ,jpcbeg,kpcbeg,iNOff,ijNOff ))
    twopi  = 2._RFREAL*global%pi

    DO k = kpcbeg-region%nDumCells, kpcend+region%nDumCells
      DO j = jpcbeg                 , jpcend
        DO i = ipcbeg-region%nDumCells, ipcend+region%nDumCells 
          ijkC        = IndIJK( i   , j   , k   ,iCOff,ijCOff) 
          iCorner(1)  = IndIJK( i   , j   , k   ,iNOff,ijNOff) 
          iCorner(2)  = IndIJK( i+1 , j   , k   ,iNOff,ijNOff) 
          iCorner(3)  = IndIJK( i   , j+1 , k   ,iNOff,ijNOff) 
          iCorner(4)  = IndIJK( i+1 , j+1 , k   ,iNOff,ijNOff) 
          iCorner(5)  = IndIJK( i   , j   , k+1 ,iNOff,ijNOff) 
          iCorner(6)  = IndIJK( i+1 , j   , k+1 ,iNOff,ijNOff) 
          iCorner(7)  = IndIJK( i   , j+1 , k+1 ,iNOff,ijNOff) 
          iCorner(8)  = IndIJK( i+1 , j+1 , k+1 ,iNOff,ijNOff) 

          xyzCell(:) = 0._RFREAL
          DO n=1,8
            xyzCell(1) = xyzCell(1) + xyz(XCOORD,iCorner(n))
            xyzCell(2) = xyzCell(2) + xyz(YCOORD,iCorner(n))
            xyzCell(3) = xyzCell(3) + xyz(ZCOORD,iCorner(n))
          ENDDO
          xyzCell(1) = 0.125_RFREAL*xyzCell(1)
          xyzCell(2) = 0.125_RFREAL*xyzCell(2)
          xyzCell(3) = 0.125_RFREAL*xyzCell(3)

! ------- X and Z momentum components

          yscale = (delta - ABS( xyzCell(2) ))/delta
          ywaves = SIN( nYwav*yscale*twopi )
          fy     = ampl*(lambda*yscale)**2/ &
                   (1._RFREAL + (lambda*yscale)**4)  ! normal envelope
          fy     = ywaves*fy                         ! normal oscillations

          cv(CV_MIXT_XMOM,ijkC) = cv(CV_MIXT_XMOM,ijkC) + &
                                  cv(CV_MIXT_DENS,ijkC)*umax* &
                                ( fy*SIN( nXwav*xyzCell(1)/extn*twopi )+ &
                                  fy*COS( nZwav*xyzCell(3)/span*twopi ) ) 
!          cv(CV_MIXT_ZMOM,ijkC) = cv(CV_MIXT_ZMOM,ijkC) + &
!                                  cv(CV_MIXT_DENS,ijkC)*umax* &
!                                ( fy*SIN( nXwav*xyzCell(1)/extn*twopi )+ &
!                                  fy*COS( nZwav*xyzCell(3)/span*twopi ) ) 

! ------- Y momentum component

          yscale = xyzCell(2)/delta
          ywaves = SIN( nYwav/2*yscale*twopi )
          yscale = (delta - ABS( xyzCell(2) ))/delta
          fy     = ampl*(lambda*yscale)**2/ &
                   (1._RFREAL + (lambda*yscale)**4)  ! normal envelope
          fy     = ywaves*fy                         ! normal oscillations

          cv(CV_MIXT_YMOM,ijkC) = cv(CV_MIXT_YMOM,ijkC) + &
                                  cv(CV_MIXT_DENS,ijkC)*umax* &
                                ( fy*COS( nXwav*xyzCell(1)/extn*twopi )+ &
                                  fy*COS( nZwav*xyzCell(3)/span*twopi ) ) 

        END DO
      END DO
    END DO 

  ELSE

! - previous pressure gradient is read from main solution file

    region%periInput%meanPgrad = global%moduleVar(1)
    write(*,*) region%procId,' channel MeanPgrad ', region%periInput%meanPgrad

  ENDIF       

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CoCnlInitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coCnlInitSolutionFlo.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/04/27 19:42:11  wasistho
! modified cnl initial sinusoidal perturbations
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.5  2004/04/01 03:38:23  wasistho
! modified channel initial perturbations
!
! Revision 1.4  2004/03/31 00:43:21  wasistho
! corrected variable extn
!
! Revision 1.3  2004/03/30 21:22:18  wasistho
! add perturbation on initial laminar channel flow
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************








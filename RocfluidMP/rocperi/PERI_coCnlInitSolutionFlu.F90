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
! $Id: PERI_coCnlInitSolutionFlu.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
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

! ... parameters
  TYPE (t_region) :: region

! ... loop variables
  INTEGER :: ijkC, ijkV

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER  :: ijkCr, ijkN, ijkN1
  REAL(RFREAL), POINTER :: xyz(:,:), cofg(:,:), cv(:,:)
  REAL(RFREAL) :: rgas, yc, muel, refMeanPgrad, delta
  REAL(RFREAL) :: rho, pres, Vm, Eo, xyzCell(3), twopi, xMin, xMax, zMin, zMax
  REAL(RFREAL) :: lambda,ampl,umax,nXwav,nYwav,nZwav,span,extn,fy,yscale,ywaves

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coCnlInitSolutionFlu.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_CoCnlInitSolution',&
  'PERI_coCnlInitSolutionFlu.F90' )

! get parameters and pointers ----------------------------------------------

  xyz  => region%grid%xyz
  cofg => region%grid%cofg
  cv   => region%mixt%cv

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

    DO ijkC = 1, region%grid%nCellsTot

      yc = cofg(YCOORD,ijkC)
      cv(CV_MIXT_XMOM,ijkC) = -0.5_RFREAL/muel*refMeanPgrad* &
                               (delta*delta - yc*yc)*rho

      Vm = SQRT( cv(CV_MIXT_XMOM,ijkC)*cv(CV_MIXT_XMOM,ijkC) + &
                 cv(CV_MIXT_YMOM,ijkC)*cv(CV_MIXT_YMOM,ijkC) + &	
                 cv(CV_MIXT_ZMOM,ijkC)*cv(CV_MIXT_ZMOM,ijkC))/rho
      Eo = MixtPerf_Eo_DGPVm(rho,global%refGamma,pres,Vm)
      cv(CV_MIXT_ENER,ijkC) = rho*Eo
    END DO

! - add perturbations, first define xMin/Max, zMin/Max

    xMin =  100000._RFREAL
    xMax = -100000._RFREAL
    zMin =  100000._RFREAL
    zMax = -100000._RFREAL
    DO ijkV = 1,region%grid%nVertTot
      xMin = MIN( xMin, xyz(XCOORD,ijkV) )
      xMax = MAX( xMax, xyz(XCOORD,ijkV) )
      zMin = MIN( zMin, xyz(ZCOORD,ijkV) )
      zMax = MAX( zMax, xyz(ZCOORD,ijkV) )
    ENDDO

    lambda = 5._RFREAL
    ampl   = 0.2_RFREAL
    umax   = -0.5_RFREAL/muel*refMeanPgrad*delta*delta
    nXwav  = 4._RFREAL
    nYwav  = 2._RFREAL
    nZwav  = 2._RFREAL
    span   = zMax-zMin
    extn   = xMax-xMin
    twopi  = 2._RFREAL*global%pi

    DO ijkC = 1,region%grid%nCellsTot

! --- X and Z momentum components

      xyzCell(1:3) = cofg(XCOORD:ZCOORD,ijkC)

      yscale = (delta - ABS( xyzCell(2) ))/delta
      ywaves = SIN( nYwav*yscale*twopi )
      fy     = ampl*(lambda*yscale)**2/ &
               (1._RFREAL + (lambda*yscale)**4)  ! normal envelope
      fy     = ywaves*fy                         ! normal oscillations

      cv(CV_MIXT_XMOM,ijkC) = cv(CV_MIXT_XMOM,ijkC) + &
                              cv(CV_MIXT_DENS,ijkC)*umax* &
                            ( fy*SIN( nXwav*xyzCell(1)/extn*twopi )+ &
                              fy*SIN( nZwav*xyzCell(3)/span*twopi ) ) 
      cv(CV_MIXT_ZMOM,ijkC) = cv(CV_MIXT_ZMOM,ijkC) + &
                              cv(CV_MIXT_DENS,ijkC)*umax* &
                            ( fy*SIN( nXwav*xyzCell(1)/extn*twopi )+ &
                              fy*SIN( nZwav*xyzCell(3)/span*twopi ) ) 

! --- Y momentum component

!      yscale = xyzCell(2)/delta
!      ywaves = SIN( nYwav/2*yscale*twopi )
!      yscale = (delta - ABS( xyzCell(2) ))/delta
!      fy     = ampl*(lambda*yscale)**2/ &
!               (1._RFREAL + (lambda*yscale)**4)  ! normal envelope
!      fy     = ywaves*fy                         ! normal oscillations

!      cv(CV_MIXT_YMOM,ijkC) = cv(CV_MIXT_YMOM,ijkC) + &
!                              cv(CV_MIXT_DENS,ijkC)*umax* &
!                            ( fy*COS( nXwav*xyzCell(1)/extn*twopi )+ &
!                              fy*COS( nZwav*xyzCell(3)/span*twopi ) ) 

    END DO

  ELSE

! - previous pressure gradient is read from main solution file

    region%periInput%meanPgrad = global%moduleVar(1)
    write(*,*) global%myProcId,' channel MeanPgrad ', region%periInput%meanPgrad

  ENDIF       

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CoCnlInitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coCnlInitSolutionFlu.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/17 20:03:14  wasistho
! compiled with RFLU
!
! Revision 1.2  2004/06/17 00:50:09  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/09 01:10:26  wasistho
! changed nomenclature
!
!
!
!******************************************************************************








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
! Purpose: Initialisation of CPR solutions
!
! Description: IF starts from t=0, set density and pressure as single values
!              as function of axial position using 1D compressible relations,
!              can be found in AIAA-86-1447 (Traineau, Hervat and Kuentsmann).
!              The remaining conservative variables are due to Taylor solution
!              which vary in injection normal direction.
!              IF restart from t/=0, read solution from main solution file
!              and cpr pressure gradient from specific cpr solution file.
!
! Input: region = data of current region
!
! Output: CPR solution to start simulation
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_coCprInitSolutionFlu.F90,v 1.5 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CoCprInitSolution( region )

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
  INTEGER :: ijkC

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER     :: global

  INTEGER  :: ijkCr, ijkN, ijkN1
  REAL(RFREAL), POINTER :: xyz(:,:), cofg(:,:), cv(:,:), dv(:,:)
  REAL(RFREAL) :: pi, gamma, gogampls, delta, minj, mfRat, refMassFlux
  REAL(RFREAL) :: throatMFlux, cprEps, headPres, headTemp, headSv, rgas
  REAL(RFREAL) :: phi, mach, yc, denom, choklen, hvinj
  REAL(RFREAL) :: rho, pres, Vm, Eo

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coCprInitSolutionFlu.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_CoCprInitSolution',&
  'PERI_coCprInitSolutionFlu.F90' )

! get parameters and pointers ----------------------------------------------

  xyz  => region%grid%xyz
  cofg => region%grid%cofg
  cv   => region%mixt%cv
  dv   => region%mixt%dv

  pi       = global%pi
  rgas     = 287._RFREAL
  gamma    = global%refGamma
  gogampls = gamma / (gamma + 1._RFREAL)
  cprEps   = region%periInput%cprEpsilon
  headPres = region%periInput%headPres
  headTemp = region%periInput%headTemp
  delta    = global%refLength
write(*,*)global%currentTime
  IF ((global%flowType==FLOW_UNSTEADY .AND. &
       global%currentTime < PERI_REAL_SMALL) .OR. &
      (global%flowType==FLOW_STEADY .AND. global%currentIter==0)) THEN

! - first cleanup conservative variables cv     
    cv(:,:) = 0._RFREAL
write(*,*)headPres,headTemp
! - refMassFlux is the mass flux at the streamwise location x/h=1/cprEpsilon
    minj        = region%periInput%minjRate 
    refMassFlux = region%periInput%bulkmFlux

! - bulk averaged mass flux corresponding to CPR sonic conditions
    throatMFlux = 2._RFREAL*delta*gamma*headPres/ &
                  SQRT(2._RFREAL*gogampls*rgas*headTemp)

! - non dimensinal distance head end to choking location L/h:
    headSv  = SQRT( gamma*rgas*headTemp )
    denom   = headSv*minj*SQRT( 2._RFREAL*(gamma+1._RFREAL) )
    choklen = headPres*gamma/denom

! - mfRat is ratio of bulk mass flux to mass flux corresp. to sonic conditions
    mfRat = refMassFlux/throatMFlux

! - or mfRat can be defined as axial coordinate normalized by choklen
    mfRat = 1._RFREAL/cprEps/chokLen

! - compute nondimensional axial parameter phi defined as
    phi = SQRT(1._RFREAL - mfRat*mfRat)

! - pressure and density axial distribution
    rho  = 0.5_RFREAL*headPres*(1._RFREAL+phi)/(rgas*headTemp)
    pres = (1.d0+gamma*phi)*headPres/(1._RFREAL+gamma)
    cv(CV_MIXT_DENS,:) = rho
    dv(DV_MIXT_PRES,:) = pres
    cv(CV_MIXT_ZMOM,:) = 0._RFREAL

! - mach number and head-end injection velocity
    mach = SQRT((1._RFREAL-phi)/(1._RFREAL+gamma*phi))
    hvinj= minj*rgas*headTemp/headPres
write(*,*)cv(CV_MIXT_DENS,1),dv(DV_MIXT_PRES,1),choklen,mach
write(*,*)mfRat,phi,mach,hvinj

    DO ijkC = 1, region%grid%nCellsTot
      yc = cofg(YCOORD,ijkC)

      cv(CV_MIXT_XMOM,ijkC) = minj/rho*0.5_RFREAL*pi/cprEps* &        
                              COS( 0.5_RFREAL*pi*yc/delta )*rho
      cv(CV_MIXT_YMOM,ijkC) = -minj/rho*SIN( 0.5_RFREAL*pi*yc/delta )*rho

      Vm = SQRT( cv(CV_MIXT_XMOM,ijkC)*cv(CV_MIXT_XMOM,ijkC) + &
                 cv(CV_MIXT_YMOM,ijkC)*cv(CV_MIXT_YMOM,ijkC) + &	
                 cv(CV_MIXT_ZMOM,ijkC)*cv(CV_MIXT_ZMOM,ijkC))/rho
      Eo = MixtPerf_Eo_DGPVm(rho,gamma,pres,Vm)
      cv(CV_MIXT_ENER,ijkC) = rho*Eo
    END DO

  ELSE

! - previous pressure gradient is read from main solution file

    region%periInput%meanPgrad = global%moduleVar(1)
    write(*,*) global%myProcId,' cprMeanPgrad ', region%periInput%meanPgrad

  ENDIF       

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CoCprInitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coCprInitSolutionFlu.F90,v $
! Revision 1.5  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/17 20:03:27  wasistho
! compiled with RFLU
!
! Revision 1.2  2004/06/17 00:50:28  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/09 01:10:26  wasistho
! changed nomenclature
!
!
!
!******************************************************************************










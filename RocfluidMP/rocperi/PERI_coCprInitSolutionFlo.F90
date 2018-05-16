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
! $Id: PERI_coCprInitSolutionFlo.F90,v 1.3 2008/12/06 08:44:36 mtcampbe Exp $
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

#include "Indexing.h"

! ... parameters
  TYPE (t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, jr

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER     :: global

  INTEGER  :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER  :: iCOff, ijCOff, iNOff, ijNOff, ijkC, ijkCr, ijkN, ijkN1
  REAL(RFREAL), POINTER :: xyz(:,:), cv(:,:), dv(:,:)
  REAL(RFREAL) :: pi, gamma, gogampls, delta, minj, mfRat, refMassFlux
  REAL(RFREAL) :: throatMFlux, cprEps, headPres, headTemp, headSv, rgas
  REAL(RFREAL) :: phi, mach, yp, ym, yc, dy, denom, choklen, hvinj
  REAL(RFREAL) :: rho, pres, Vm, Eo

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coCprInitSolutionFlo.F90,v $ $Revision: 1.3 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_CoCprInitSolution',&
  'PERI_coCprInitSolutionFlo.F90' )

! get parameters and pointers ----------------------------------------------

  iLev   =  region%currLevel
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz    => region%levels(iLev)%grid%xyz
  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv

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
    dv(DV_MIXT_WVEL,:) = 0._RFREAL

! - mach number and head-end injection velocity
    mach = SQRT((1._RFREAL-phi)/(1._RFREAL+gamma*phi))
    hvinj= minj*rgas*headTemp/headPres
write(*,*)cv(CV_MIXT_DENS,1),dv(DV_MIXT_PRES,1),choklen,mach
write(*,*)mfRat,phi,mach,hvinj

    DO k = kpcbeg, kpcbeg
      DO i = ipcbeg, ipcbeg
        DO j = jpcbeg-region%nDumCells, jpcend+region%nDumCells
          ijkC  = IndIJK( i , j  , k ,iCOff,ijCOff) 
          jr    = MAX( j, jpcbeg )
          jr    = MIN( jr,jpcend )
          ijkN  = IndIJK( i , jr  , k ,iNOff,ijNOff) 
          ijkN1 = IndIJK( i , jr+1, k ,iNOff,ijNOff) 
          dy    = xyz(YCOORD,ijkN1)-xyz(YCOORD,ijkN)
          yp    = xyz(YCOORD,ijkN1)
          ym    = xyz(YCOORD,ijkN) 
          yc    = 0.5_RFREAL*(xyz(YCOORD,ijkN1)+xyz(YCOORD,ijkN))
          IF (dy < 0._RFREAL) THEN
            CALL ErrorStop( region%global,ERR_PERI_GEO,__LINE__, &
                            'grid contains negative spacing' )
          ENDIF
            
!          dv(DV_MIXT_UVEL,ijkC) = 0.5_RFREAL*pi*refMassFlux* &
!                                  (SIN(0.5_RFREAL*pi*yp/delta)- &
!                                  SIN(0.5_RFREAL*pi*ym/delta))/(pi*dy*rho)
!          dv(DV_MIXT_VVEL,ijkC) = cprEps*refMassFlux*(COS(0.5_RFREAL*pi*yp/ &
!                                  delta)-COS(0.5_RFREAL*pi*ym/delta))/(pi*dy*rho)

          dv(DV_MIXT_UVEL,ijkC) = minj/rho*0.5_RFREAL*pi/cprEps* &        
                                  COS( 0.5_RFREAL*pi*yc/delta )
          dv(DV_MIXT_VVEL,ijkC) = -minj/rho*SIN( 0.5_RFREAL*pi*yc/delta )

          Vm = SQRT( dv(DV_MIXT_UVEL,ijkC)*dv(DV_MIXT_UVEL,ijkC) + &
                     dv(DV_MIXT_VVEL,ijkC)*dv(DV_MIXT_VVEL,ijkC) + &	
                     dv(DV_MIXT_WVEL,ijkC)*dv(DV_MIXT_WVEL,ijkC))
          Eo = MixtPerf_Eo_DGPVm(rho,gamma,pres,Vm)
          cv(CV_MIXT_XMOM,ijkC) = rho*dv(DV_MIXT_UVEL,ijkC)
          cv(CV_MIXT_YMOM,ijkC) = rho*dv(DV_MIXT_VVEL,ijkC)
          cv(CV_MIXT_ZMOM,ijkC) = rho*dv(DV_MIXT_WVEL,ijkC)
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
!IF (i==ipcbeg.AND.k==kpcbeg) &
!write(*,*)j,cv(2,ijkC),cv(3,ijkC),cv(4,ijkC),cv(5,ijkC)
IF(j==1.AND.k==2) write(*,*)i,j,k,cv(1:5,ijkC)
        END DO
      END DO
    END DO 

  ELSE

! - previous pressure gradient is read from main solution file

    region%periInput%meanPgrad = global%moduleVar(1)
    write(*,*) region%procId,' cprMeanPgrad ', region%periInput%meanPgrad

  ENDIF       

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CoCprInitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coCprInitSolutionFlo.F90,v $
! Revision 1.3  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
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








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
! ******************************************************************************
!
! Purpose: Finish computation of substantial derivative.
!
! Description: None.
!
! Input:
!   region             Data of current region
!
! Output: None.
!
! Notes:
!   1. On entry: region%mixt%sd     Residual from convective terms only
!                region%mixt%rhs    Residual from all terms
!   2. On exit:  region%mixt%sd     D{u,v,w}/Dt
!   3. Current implementation assumes grid is not moving
!
! ******************************************************************************
!
! $Id: RFLU_FinishSD.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_FinishSD(region)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global

! TEMPORARY
!  USE RFLU_ModTime, ONLY: RFLU_GetTimeRK
! END TEMPORARY

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: arrayLimLow,arrayLimUpp,icg
  REAL(RFREAL) :: term 
  REAL(RFREAL), DIMENSION(:), POINTER :: vol
  REAL(RFREAL), DIMENSION(:,:), POINTER :: sd,rhs,cv,tv
  TYPE(t_global), POINTER :: global

! TEMPORARY
!  REAL(RFREAL) :: term2,term3,DuDt,DuDt2,DuDt3,xr,xl,dx
! END TEMPORARY

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'RFLU_FinishSD',&
  'RFLU_FinishSD.F90')

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************

  IF ( region%mixtInput%moveGrid .EQV. .TRUE. ) THEN
    CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
      'Equilibrium Eulerian method not yet implemented with moving grid')
  END IF ! region%mixtInput%moveGrid

  IF ( region%mixtInput%indSd /= 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! region%mixtInput%indSd

! ******************************************************************************
! Set up pointers and compute array bounds
! ******************************************************************************

  vol => region%grid%vol

  sd  => region%mixt%sd
  rhs => region%mixt%rhs
  cv  => region%mixt%cv
  tv  => region%mixt%tv

! ******************************************************************************
! Compute DuDt 
! ******************************************************************************

! DEBUG
!  dx = 2.0_RFREAL/(region%grid%nCells/9.0_RFREAL)
!  WRITE(*,*) 'RFLU_FinishSD:',dx
! END DEBUG

  DO icg = 1,region%grid%nCellsTot
    term = 1.0_RFREAL/(cv(CV_MIXT_DENS,icg)*vol(icg))    

    sd(SD_XMOM,icg) = term*(sd(SD_XMOM,icg) - rhs(CV_MIXT_XMOM,icg))
    sd(SD_YMOM,icg) = term*(sd(SD_YMOM,icg) - rhs(CV_MIXT_YMOM,icg))
    sd(SD_ZMOM,icg) = term*(sd(SD_ZMOM,icg) - rhs(CV_MIXT_ZMOM,icg))

! TEMPORARY - Comparison with shock tube Sod 1 at 1ms
!    term = region%grid%cofg(XCOORD,icg)/1.0E-3_RFREAL 
!        
!    IF ( (ABS(region%grid%cofg(YCOORD,icg)) < 0.017_RFREAL) .AND. & 
!         (ABS(region%grid%cofg(ZCOORD,icg)) < 0.017_RFREAL) ) THEN  
!    IF ( region%grid%cofg(XCOORD,icg) > -0.374165735491392_RFREAL .AND. & 
!         region%grid%cofg(XCOORD,icg) < -2.222221147508515E-002 ) THEN 
!      term2 = -5.0_RFREAL/6.0_RFREAL*term/1.0E-3_RFREAL             
!      term3 = 25.0_RFREAL/36.0_RFREAL*(374.166_RFREAL+term)/1.0E-3_RFREAL   
! 
!      xl = region%grid%cofg(XCOORD,icg)-0.5_RFREAL*dx
!      xr = region%grid%cofg(XCOORD,icg)+0.5_RFREAL*dx
!      DuDt = ((311.805_RFREAL-xr/6.0E-3_RFREAL)**7-(311.805_RFREAL-xl/6.0E-3_RFREAL)**7)/ &
!             ((311.805_RFREAL-xr/6.0E-3_RFREAL)**6-(311.805_RFREAL-xl/6.0E-3_RFREAL)**6)/(1.4_RFREAL*1.0E-3_RFREAL)
!      DuDt2 = 5.0_RFREAL/6.0_RFREAL*(311.805_RFREAL-region%grid%cofg(XCOORD,icg)/(6.0_RFREAL*1.0E-3_RFREAL))/1.0E-3_RFREAL  
!      DuDt3 = 5.0_RFREAL/(6.0_RFREAL*dx)*((xr/1.0E-3_RFREAL*(311.805-xr/(12.0_RFREAL*1.0E-3_RFREAL))) & 
!                                        - (xl/1.0E-3_RFREAL*(311.805-xl/(12.0_RFREAL*1.0E-3_RFREAL))))
!    ELSE 
!      term2 = 0.0_RFREAL
!      term3 = 0.0_RFREAL
!      DuDt  = 0.0_RFREAL
!      DuDt2 = 0.0_RFREAL
!      DuDt3 = 0.0_RFREAL
!    END IF 
!
!    WRITE(11,'(7(1X,E23.16))') region%grid%cofg(XCOORD,icg), & 
!                               sd(SD_XMOM,icg), &
!                               term2+term3, &
!                               DuDt,DuDt2,DuDt3, & 
!                               region%mixt%cv(CV_MIXT_XMOM,icg)/region%mixt%cv(CV_MIXT_DENS,icg)
!    END IF 
! END TEMPORARY                   
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_FinishSD

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_FinishSD.F90,v $
! Revision 1.4  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/03/31 16:49:32  haselbac
! Improved computation of sd
!
! Revision 1.1  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! ******************************************************************************








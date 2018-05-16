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
! Purpose: compute dependent variables for a thermally and calorically
!          perfect gas.
!
! Description: none.
!
! Input: inBeg  = first value to update
!        inEnd  = last value to update
!        indCp  = indicates if cp varies over cells (=1) or is constant (=0)
!        indMol = indicates if the mol mass varies over cells (=1) or is
!                 constant (=0)
!        cv     = conservative variables
!        gv     = gas variables (cp, Mol)
!
! Output: dv = dependent variables (p, T, c and u, v, w for RocfloMP)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PerfgasDependentVars.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PerfgasDependentVars( inBeg,inEnd,indCp,indMol,cv,gv,dv )

  USE ModDataTypes
  USE ModParameters

#ifdef RFLU
  USE ModInterfaces, ONLY: MixtPerf_C_GRT,MixtPerf_G_CpR, & 
                           MixtPerf_P_DEoGVm2,MixtPerf_R_M,MixtPerf_T_DPR
#endif

  IMPLICIT NONE

! ... parameters
  INTEGER :: inBeg, inEnd, indCp, indMol

  REAL(RFREAL), POINTER :: cv(:,:), gv(:,:), dv(:,:)

! ... loop variables
  INTEGER :: ic

! ... local variables
  REAL(RFREAL) :: rgas,rrho,Vm2
#ifdef RFLO
  REAL(RFREAL) :: gam1, g1cp
#endif
#ifdef RFLU
  REAL(RFREAL) :: Eo,gamma,rho
#endif

!******************************************************************************

  DO ic=inBeg,inEnd
#ifdef RFLO
    rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ic*indMol)
    gam1 = gv(GV_MIXT_CP,ic*indCp)/(gv(GV_MIXT_CP,ic*indCp)-rgas) - 1._RFREAL
    g1cp = gam1*gv(GV_MIXT_CP,ic*indCp)
    rrho = 1._RFREAL/cv(CV_MIXT_DENS,ic)

    dv(DV_MIXT_UVEL,ic) = cv(CV_MIXT_XMOM,ic)*rrho
    dv(DV_MIXT_VVEL,ic) = cv(CV_MIXT_YMOM,ic)*rrho
    dv(DV_MIXT_WVEL,ic) = cv(CV_MIXT_ZMOM,ic)*rrho

    vm2 = dv(DV_MIXT_UVEL,ic)*dv(DV_MIXT_UVEL,ic) + &
          dv(DV_MIXT_VVEL,ic)*dv(DV_MIXT_VVEL,ic) + &
          dv(DV_MIXT_WVEL,ic)*dv(DV_MIXT_WVEL,ic)

    dv(DV_MIXT_PRES,ic) = gam1*(cv(CV_MIXT_ENER,ic)- &
                                0.5_RFREAL*vm2*cv(CV_MIXT_DENS,ic))
    dv(DV_MIXT_TEMP,ic) = dv(DV_MIXT_PRES,ic)*rrho/rgas
    dv(DV_MIXT_SOUN,ic) = SQRT(g1cp*dv(DV_MIXT_TEMP,ic))
#endif
#ifdef RFLU
    rgas  = MixtPerf_R_M(gv(GV_MIXT_MOL,ic*indMol))
    gamma = MixtPerf_G_CpR(gv(GV_MIXT_CP,ic*indCp),rgas)

    rho  = cv(CV_MIXT_DENS,ic)
    rrho = 1.0_RFREAL/rho
    Eo   = cv(CV_MIXT_ENER,ic)*rrho

    Vm2 = (cv(CV_MIXT_XMOM,ic)*cv(CV_MIXT_XMOM,ic) + &
           cv(CV_MIXT_YMOM,ic)*cv(CV_MIXT_YMOM,ic) + &
           cv(CV_MIXT_ZMOM,ic)*cv(CV_MIXT_ZMOM,ic))*rrho*rrho

    dv(DV_MIXT_PRES,ic) = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
    dv(DV_MIXT_TEMP,ic) = MixtPerf_T_DPR(rho,dv(DV_MIXT_PRES,ic),rgas)
    dv(DV_MIXT_SOUN,ic) = MixtPerf_C_GRT(gamma,rgas,dv(DV_MIXT_TEMP,ic))
#endif
  ENDDO

END SUBROUTINE PerfgasDependentVars

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PerfgasDependentVars.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:49:56  haselbac
! Initial revision after changing case
!
! Revision 1.12  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/05/13 23:46:52  haselbac
! Reverted to use of MixtPerf routines for RFLU
!
! Revision 1.8  2003/05/06 20:05:39  jblazek
! Corrected bug in grid motion (corner "averaging").
!
! Revision 1.7  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.6  2002/06/05 18:33:00  haselbac
! Converted to use of mixtPerf routines
!
! Revision 1.5  2002/03/18 22:25:45  jblazek
! Finished multiblock and MPI.
!
! Revision 1.4  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.2  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/10 00:02:06  jblazek
! Added calculation of mixture properties.
!
!******************************************************************************







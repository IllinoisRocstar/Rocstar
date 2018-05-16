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
! Purpose: compute additional terms to NS due to CPR formulations
!
! Description: all variables in slow terms are averaged over i,k plane.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = CPR slow terms added to the residual.
!
! Notes: This routine contents MPI global SUM which does not account for
!        summation over regions in same processor. Therefore nProcAlloc should
!        be equal to nRegions.
!
!******************************************************************************
!
! $Id: PERI_coCprSlowTermsFlo.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CoCprSlowTerms( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            MixtPerf_R_M, MixtPerf_G_CpR
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, ijkC0, nVar, ipc, jpc, kpc
  INTEGER :: indCp, indMol

  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), rhs(:,:), vol(:)
  REAL(RFREAL), POINTER :: cprVar(:,:), varSend(:,:), varRecv(:,:)
  REAL(RFREAL), ALLOCATABLE :: rgas(:), cpgas(:)
  REAL(RFREAL) :: volEps, averSize, sndAverSize, rcvAverSize, rAvgRho
  REAL(RFREAL) :: denumer, denomin, delta, rDelta, gamma, gammin, gagmin
  REAL(RFREAL) :: rho, ruc, ren, uve, vve, wve, tmp, prs, rPrim, uPrim, vPrim 
  REAL(RFREAL) :: wPrim, tPrim, pPrim, uGrad, vGrad, wGrad, tGrad, pGrad
  REAL(RFREAL) :: refMeanPgrad, fCont, fMomX, fMomY, fMomZ, fRhoE

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coCprSlowTermsFlo.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( region%global,'PERI_CoCprSlowTerms',&
  'PERI_coCprSlowTermsFlo.F90' )

! get dimensions, pointers and parameters -------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  gv     => region%levels(iLev)%mixt%gv
  rhs    => region%levels(iLev)%mixt%rhs
  vol    => region%levels(iLev)%grid%vol
  cprVar => region%levels(iLev)%peri%cprVar
  varSend=> region%levels(iLev)%peri%varSend
  varRecv=> region%levels(iLev)%peri%varRecv

  ipc      = ipcend-ipcbeg+1
  jpc      = jpcend-jpcbeg+1
  kpc      = kpcend-kpcbeg+1
  averSize = DBLE(ipc*kpc)

  indCp        = region%levels(iLev)%mixt%indCp
  indMol       = region%levels(iLev)%mixt%indMol
  nVar         = region%periInput%nVar

  refMeanPgrad = region%periInput%meanPgrad
  delta        = region%global%refLength
  rDelta       = 1._RFREAL/delta

! allocate temporary space

  ALLOCATE( rgas(jpcbeg:jpcend), cpgas(jpcbeg:jpcend) )

! x-z plane averaged variables

  DO j=jpcbeg,jpcend
    cprVar(1:nVar,j) = 0._RFREAL
    rgas(j)          = 0._RFREAL
    cpgas(j)         = 0._RFREAL
    DO k=kpcbeg,kpcend
      DO i=ipcbeg,ipcend
        ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)

        cprVar(CPR_RHO,j) = cprVar(CPR_RHO,j) + cv(CV_MIXT_DENS,ijkC0)
        cprVar(CPR_RUC,j) = cprVar(CPR_RUC,j) + cv(CV_MIXT_XMOM,ijkC0)
        cprVar(CPR_RVC,j) = cprVar(CPR_RVC,j) + cv(CV_MIXT_YMOM,ijkC0)
        cprVar(CPR_UVE,j) = cprVar(CPR_UVE,j) + dv(DV_MIXT_UVEL,ijkC0)
        cprVar(CPR_VVE,j) = cprVar(CPR_VVE,j) + dv(DV_MIXT_VVEL,ijkC0)
        cprVar(CPR_TMP,j) = cprVar(CPR_TMP,j) + dv(DV_MIXT_TEMP,ijkC0)
        cprVar(CPR_PRS,j) = cprVar(CPR_PRS,j) + dv(DV_MIXT_PRES,ijkC0)

        rgas(j) = rgas(j) + MixtPerf_R_M( gv(GV_MIXT_MOL,ijkC0*indMol) )
        cpgas(j)= cpgas(j)+ gv(GV_MIXT_CP,ijkC0*indCp)
      ENDDO
    ENDDO
  ENDDO

#ifdef MPI
  DO j = jpcbeg, jpcend
    DO l = 1,nVar
      varSend(j-jpcbeg+1,l) = cprVar(l,j)
    END DO
    varSend(j-jpcbeg+1,nVar+1) = rgas(j)
    varSend(j-jpcbeg+1,nVar+2) = cpgas(j)
  END DO
  CALL MPI_BARRIER( region%global%mpiComm, region%global%mpierr )
  DO l = 1, nVar+GAS_NVAR
    CALL MPI_ALLREDUCE( varSend(1,l),varRecv(1,l),jpc,MPI_DOUBLE_PRECISION, &
         MPI_SUM,region%global%mpiComm, region%global%mpierr )
  END DO
  DO j = jpcbeg, jpcend
    DO l = 1,nVar
      cprVar(l,j) = varRecv(j-jpcbeg+1,l)
    END DO
    rgas(j)  = varRecv(j-jpcbeg+1,nVar+1)
    cpgas(j) = varRecv(j-jpcbeg+1,nVar+2)
  END DO

  sndAverSize = averSize
  CALL MPI_ALLREDUCE( sndAverSize,rcvAverSize,1, MPI_DOUBLE_PRECISION, &
       MPI_SUM,region%global%mpiComm,region%global%mpierr )
  averSize = rcvAverSize
#endif

  DO j = jpcbeg,jpcend
    DO l = 1,nVar
      cprVar(l,j) = cprVar(l,j)/averSize
    ENDDO
    rgas(j) = rgas(j)/averSize
    cpgas(j)= cpgas(j)/averSize

    gamma   = MixtPerf_G_CpR( cpgas(j),rgas(j) )
    gagmin  = gamma/(gamma - 1_RFREAL)

    rAvgRho = 1._RFREAL/cprVar(CPR_RHO,j)  
    denumer = gagmin*refMeanPgrad + (cprVar(CPR_RUC,j)**2 + &
                                     cprVar(CPR_RVC,j)**2)*rDelta*rAvgRho
    denomin = gagmin*cprVar(CPR_PRS,j)*rAvgRho + &
              (cprVar(CPR_RUC,j)**2+cprVar(CPR_RVC,j)**2)*rAvgRho**2
    cprVar(CPR_DOR,j) = denumer/denomin
    cprVar(CPR_DOU,j) = (cprVar(CPR_RUC,j)*rDelta-cprVar(CPR_UVE,j)* &
                        cprVar(CPR_DOR,j))*rAvgRho
    cprVar(CPR_DOV,j) = -cprVar(CPR_VVE,j)*rAvgRho*cprVar(CPR_DOR,j)
    cprVar(CPR_DOT,j) = (refMeanPgrad/rgas(j)-cprVar(CPR_TMP,j)* &
                        cprVar(CPR_DOR,j))*rAvgRho
  ENDDO

  DO k = kpcbeg, kpcend
    DO j = jpcbeg, jpcend
      DO i = ipcbeg, ipcend
        ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)

        rho = cv(CV_MIXT_DENS,ijkC0)
        ruc = cv(CV_MIXT_XMOM,ijkC0)
        ren = cv(CV_MIXT_ENER,ijkC0)
        uve = dv(DV_MIXT_UVEL,ijkC0)
        vve = dv(DV_MIXT_VVEL,ijkC0)
        wve = dv(DV_MIXT_WVEL,ijkC0)
        tmp = dv(DV_MIXT_TEMP,ijkC0)
        prs = dv(DV_MIXT_PRES,ijkC0)

        rPrim = rho - cprVar(CPR_RHO,j)
        uPrim = uve - cprVar(CPR_UVE,j)
        vPrim = vve - cprVar(CPR_VVE,j)
        wPrim = wve
        tPrim = tmp - cprVar(CPR_TMP,j)
        pPrim = prs - cprVar(CPR_PRS,j)

        uGrad = uPrim*rDelta + cprVar(CPR_DOU,j)
        vGrad = vPrim*rDelta + cprVar(CPR_DOV,j)  
        wGrad = wPrim*rDelta
        tGrad = tPrim*rDelta + cprVar(CPR_DOT,j)  
        pGrad = pPrim*rDelta + refMeanPgrad 

        fCont = (cprVar(CPR_RUC,j)+cprVar(CPR_RHO,j)*uPrim+cprVar(CPR_UVE,j)* &
                rPrim)*rDelta + uPrim*cprVar(CPR_DOR,j) + rPrim*cprVar(CPR_DOU,j)
        fMomX = uve*fCont + ruc*uGrad + pGrad
        fMomY = vve*fCont + ruc*vGrad
        fMomZ = wve*fCont + ruc*wGrad
        gamma = MixtPerf_G_CpR( cpgas(j),rgas(j) )
        gammin= gamma-1._RFREAL
        fRhoE = ren/rho*fCont + ruc*(uve*uGrad+vve*vGrad+wve*wGrad+rgas(j)/ &
                gammin*tGrad) + prs*uGrad +uve*pGrad

        volEps= vol(ijkC0)*region%periInput%cprEpsilon
        rhs(CV_MIXT_DENS,ijkC0) = rhs(CV_MIXT_DENS,ijkC0) + volEps*fCont
        rhs(CV_MIXT_XMOM,ijkC0) = rhs(CV_MIXT_XMOM,ijkC0) + volEps*fMomX
        rhs(CV_MIXT_YMOM,ijkC0) = rhs(CV_MIXT_YMOM,ijkC0) + volEps*fMomY
        rhs(CV_MIXT_ZMOM,ijkC0) = rhs(CV_MIXT_ZMOM,ijkC0) + volEps*fMomZ
        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + volEps*fRhoE
!=======================================
!        fCont = cprVar(CPR_RUC,j)*rdelta
!        fMomX = cprVar(CPR_UVE,j)*fCont 
!        fMomY = cprVar(CPR_VVE,j)*fCont
!        fMomZ = 0._RFREAL
!        fRhoE = ren/rho*fCont + ruc*(uve*uGrad+vve*vGrad+wve*wGrad+rgas(j)/ &
!                gammin*tGrad) + prs*uGrad +uve*pGrad
!        rhs(CV_MIXT_DENS,ijkC0) = rhs(CV_MIXT_DENS,ijkC0) + volEps*fCont
!        rhs(CV_MIXT_XMOM,ijkC0) = rhs(CV_MIXT_XMOM,ijkC0) + volEps* &
!                                  (fMomX + ruc*uve*rdelta + refMeanPgrad)
!        rhs(CV_MIXT_YMOM,ijkC0) = rhs(CV_MIXT_YMOM,ijkC0) + volEps*ruc*vve*rdelta
!        rhs(CV_MIXT_ZMOM,ijkC0) = rhs(CV_MIXT_ZMOM,ijkC0) + volEps*ruc*wve*rdelta
!        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + volEps* &
!                                  (ren/rho+prs)*uve*rdelta
!========================================
      ENDDO
    ENDDO
  ENDDO

! deallocate temporary workspace

  DEALLOCATE( rgas, cpgas )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE PERI_CoCprSlowTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coCprSlowTermsFlo.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
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
!******************************************************************************








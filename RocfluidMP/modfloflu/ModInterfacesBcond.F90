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
! Purpose: set explicit interfaces to subroutines and functions
!          related to boundary conditions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesBcond.F90,v 1.5 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesBcond

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE BcondFarfieldPerf( machInf,alphaInf,betaInf,pInf,tInf, & 
                                sxn,syn,szn,cpgas,mol, & 
                                rho,rhou,rhov,rhow,rhoe,press, &
                                rhob,rhoub,rhovb,rhowb,rhoeb,pb )
    USE ModDataTypes
    REAL(RFREAL) :: machInf, alphaInf, betaInf, pInf, tInf
    REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
    REAL(RFREAL) :: sxn, syn, szn, cpgas, mol
    REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb, pb
  END SUBROUTINE BcondFarfieldPerf

  SUBROUTINE BcondInflowPerf( bcOptType,bcOptFixed,ptot,ttot,betah,betav, & 
                              mach,sxn,syn,szn,cpgas,mm,rl,rul,rvl,rwl,rr, &
                              rur,rvr,rwr,rer,pr )
    USE ModDataTypes
    INTEGER, INTENT(IN) :: bcOptFixed,bcOptType
    REAL(RFREAL), INTENT(IN)  :: betah, betav, cpgas, mach, mm, sxn, syn, szn, &
                                 ptot, rl, rul, rvl, rwl, ttot
    REAL(RFREAL), INTENT(OUT) :: rer, rr, rur, rvr, rwr, pr
  END SUBROUTINE BcondInflowPerf

  SUBROUTINE BcondInflowPerf_GL(bcOptType,ro,po,to,Bp,Bt,cvl,cvv,cvg,Rg,Rv,ur, &
                                vr,wr,vfgr,vfvr,vflr,temp,press,nx,ny,nz,rl, &
                                rul,rvl,rwl,rel,rgpgl,rvpvl,pl,rr,rur,rvr, &
                                rwr,rer,rgpgr,rvpvr,pr)
    USE ModDataTypes
    INTEGER, INTENT(IN) :: bcOptType
    REAL(RFREAL), INTENT(IN) :: Bp,Bt,cvg,cvl,cvv,nx,ny,nz,pl,po,press,rel, &
                                Rg,rgpgl,rl,ro,rul,Rv,rvl,rvpvl,rwl,temp,to, &
                                ur,vfgr,vflr,vfvr,vr,wr 
    REAL(RFREAL), INTENT(OUT):: pr,rer,rgpgr,rr,rur,rvr,rvpvr,rwr 
  END SUBROUTINE BcondInflowPerf_GL

  SUBROUTINE BcondInjectionPerf( distrib,minj,tinj,rhoVrel,sxn,syn,szn, &
                                 cpgas,mm,p,rhob,rhoub,rhovb,rhowb,rhoeb,pb, &
                                 uinj,vinj,winj )
    USE ModDataTypes
    INTEGER,      INTENT(IN)  :: distrib
    REAL(RFREAL), INTENT(IN)  :: cpgas, minj, mm, sxn, syn, szn, tinj, &
                                 rhoVrel(3), p
    REAL(RFREAL), INTENT(OUT) :: rhob, rhoub, rhovb, rhowb, rhoeb, pb, &
                                 uinj, vinj, winj
  END SUBROUTINE BcondInjectionPerf

  SUBROUTINE BcondOutflowPerf( bcOpt,pout,sxn,syn,szn,cpgas,mol, & 
                               rho,rhou,rhov,rhow,rhoe,press, &
                               rhob,rhoub,rhovb,rhowb,rhoeb )
    USE ModDataTypes
    INTEGER :: bcOpt
    REAL(RFREAL), INTENT(IN) :: cpgas,mol,pout,press,rho,rhou, &
                                rhov,rhow,rhoe,sxn,syn,szn 
    REAL(RFREAL), INTENT(OUT) :: rhob,rhoub,rhovb,rhowb,rhoeb
  END SUBROUTINE BcondOutflowPerf

  SUBROUTINE BcondOutflowPerf_GL(bcOpt,ro,Po,To,betaP,betaT,cvl,cvv,cvg,Rg,Rv, &
                                 pout,sxn,syn,szn,rho,rhou,rhov,rhow,rhoe, &
                                 rhogpg,rhovpv,pin,rhob,rhoub,rhovb,rhowb, &
                                 rhoeb,rhogpgb,rhovpvb)
    USE ModDataTypes
    INTEGER :: bcOpt
    REAL(RFREAL), INTENT(IN) :: betaP,betaT,cvg,cvl,cvv,pin,Po,pout,Rg,rho, &
                                rhoe,rhogpg,rhou,rhov,rhovpv,rhow,ro,Rv,sxn, &
                                syn,szn,To
    REAL(RFREAL), INTENT(OUT) :: rhob,rhoeb,rhogpgb,rhoub,rhovb,rhovpvb,rhowb
  END SUBROUTINE BcondOutflowPerf_GL

  SUBROUTINE UpdateBoundaryConditionsMP( regions, istage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: istage
  END SUBROUTINE UpdateBoundaryConditionsMP

  SUBROUTINE UpdateTbc( region,t,dt,final )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    REAL(RFREAL),   INTENT(IN)    :: t, dt
    LOGICAL,        INTENT(IN)    :: final
  END SUBROUTINE UpdateTbc

  SUBROUTINE UpdateTbcPiecewise( global,tbc,t )
    USE ModDataTypes
    USE ModBndPatch, ONLY : t_tbcvalues
    USE ModGlobal,   ONLY : t_global
    TYPE(t_global),    POINTER       :: global
    TYPE(t_tbcvalues), INTENT(INOUT) :: tbc
    REAL(RFREAL),      INTENT(IN)    :: t
  END SUBROUTINE UpdateTbcPiecewise

  SUBROUTINE UpdateTbcSinusoidal( global,tbc,t )
    USE ModDataTypes
    USE ModBndPatch, ONLY : t_tbcvalues
    USE ModGlobal,   ONLY : t_global
    TYPE(t_global),    POINTER       :: global
    TYPE(t_tbcvalues), INTENT(INOUT) :: tbc
    REAL(RFREAL),      INTENT(IN)    :: t
  END SUBROUTINE UpdateTbcSinusoidal

  SUBROUTINE UpdateTbcStochastic( region,tbc,dt )
    USE ModDataTypes
    USE ModBndPatch,   ONLY : t_tbcvalues
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region),    INTENT(INOUT) :: region
    TYPE(t_tbcvalues), INTENT(INOUT) :: tbc
    REAL(RFREAL),      INTENT(IN)    :: dt
  END SUBROUTINE UpdateTbcStochastic

  SUBROUTINE UpdateTbcWhitenoise( region,tbc )
    USE ModBndPatch,   ONLY : t_tbcvalues
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region),    INTENT(INOUT) :: region
    TYPE(t_tbcvalues), INTENT(INOUT) :: tbc
  END SUBROUTINE UpdateTbcWhitenoise

  SUBROUTINE ZeroDummyCellsMP( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE ZeroDummyCellsMP

  END INTERFACE

END MODULE ModInterfacesBcond

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesBcond.F90,v $
! Revision 1.5  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/26 20:21:52  haselbac
! Added ifs for new routines
!
! Revision 1.2  2004/01/29 22:57:24  haselbac
! Changed interface for bcondInflowPerf.F90
!
! Revision 1.1  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
!******************************************************************************







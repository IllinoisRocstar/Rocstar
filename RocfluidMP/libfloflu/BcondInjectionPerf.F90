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
! Purpose: set injection boundary condition for one cell.
!
! Description: none.
!
! Input: distrib = 0: BC values constant, 1: BC values varying
!        minj    = injection mass flux
!        tinj    = injection temperature
!        rhoVrel = density * relative velocity (V-inj - boundary movement)
!        p       = static pressure (boundary cell)
!        sx/y/zn = components of ortho-normalized face vector (outward facing)
!        cpgas   = specific heat at constant pressure (boundary cell)
!        mm      = molecular mass at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!         pb        = pressure at boundary
!         u/v/winj  = components of injection velocity
!
! Notes: this condition is valid only for thermally and calorically
!        perfect gas. It is also assumed that for distrib==BCDAT_CONSTANT the
!        relative velocity is identical to injection velocity.
!
!******************************************************************************
!
! $Id: BcondInjectionPerf.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BcondInjectionPerf( distrib,minj,tinj,rhoVrel,sxn,syn,szn, &
                               cpgas,mm,p,rhob,rhoub,rhovb,rhowb,rhoeb,pb, &
                               uinj,vinj,winj )

  USE ModDataTypes
  USE ModParameters
  USE ModInterfaces, ONLY: MixtPerf_D_PRT, MixtPerf_Eo_DGPUVW, MixtPerf_G_CpR, & 
                           MixtPerf_R_M

  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: distrib

  REAL(RFREAL), INTENT(IN)  :: cpgas, minj, tinj, rhoVrel(3), mm, &
                               sxn, syn, szn, p
  REAL(RFREAL), INTENT(OUT) :: rhob, rhoub, rhovb, rhowb, rhoeb, pb, &
                               uinj, vinj, winj

! ... local variables
  REAL(RFREAL) :: g, qb, rgas, ub, vb, wb

!******************************************************************************
! gas properties

  rgas = MixtPerf_R_M( mm )
  g    = MixtPerf_G_CpR( cpgas,rgas )

! set boundary values

  pb   = p
  rhob = MixtPerf_D_PRT( pb,rgas,tinj )

  qb   = -minj/rhob
  uinj = qb*sxn
  vinj = qb*syn
  winj = qb*szn

  IF (distrib == BCDAT_CONSTANT) THEN
    ub = uinj
    vb = vinj
    wb = winj
  ELSE
    ub = rhoVrel(1)/rhob
    vb = rhoVrel(2)/rhob
    wb = rhoVrel(3)/rhob
  ENDIF

  rhoeb = rhob * MixtPerf_Eo_DGPUVW( rhob,g,pb,ub,vb,wb )
  rhoub = rhob*ub
  rhovb = rhob*vb
  rhowb = rhob*wb

END SUBROUTINE BcondInjectionPerf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondInjectionPerf.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:00  haselbac
! Initial revision after changing case
!
! Revision 1.6  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2002/09/27 00:26:17  jblazek
! Changed for moving boundaries.
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************







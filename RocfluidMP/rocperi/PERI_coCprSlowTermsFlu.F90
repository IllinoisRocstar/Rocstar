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
! $Id: PERI_coCprSlowTermsFlu.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CoCprSlowTerms( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : MixtPerf_R_M, MixtPerf_G_CpR
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coCprSlowTermsFlu.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction( region%global,'PERI_CoCprSlowTerms',&
  'PERI_coCprSlowTermsFlu.F90' )

! under construction ----------------------------------------------------------

  CALL ErrorStop(region%global,ERR_PERI_INPUT,__LINE__, &
       'CPR slow terms routine for Rocflu is not built yet')

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE PERI_CoCprSlowTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coCprSlowTermsFlu.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/06/17 20:03:40  wasistho
! compiled with RFLU
!
! Revision 1.1  2004/06/09 01:10:26  wasistho
! changed nomenclature
!
!
!
!******************************************************************************








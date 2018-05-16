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
! Purpose: define various parameters pertinent to RfloPrep
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: PREP_ModParameters.F90,v 1.6 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE PREP_ModParameters

  USE ModDataTypes
  IMPLICIT NONE

! PREP integer parameters -----------------------------------------------------

!******************************************************************************
! inflow profiles:  
!******************************************************************************

  INTEGER, PARAMETER :: INFLO_TAYLOR_CYL  = 1, &
                        INFLO_TAYLOR_PLAN = 2, &
                        INFLO_BLLAM_CYL   = 3, &
                        INFLO_BLLAM_PLAN  = 4, &
                        INFLO_BLTURB_CYL  = 5, &
                        INFLO_BLTURB_PLAN = 6   

!******************************************************************************
! hardcoded initial conditions  
!******************************************************************************

  INTEGER, PARAMETER :: INITFLO_UNIFORM       = 0, &
                        INITFLO_PISTON_EXPAN  = 1

! PREP real parameters --------------------------------------------------------

  REAL(RFREAL), PARAMETER :: REAL_SMALL   = 1.E-16_RFREAL, &
                             REAL_LARGE   = 1.E+32_RFREAL

END MODULE PREP_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_ModParameters.F90,v $
! Revision 1.6  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/09/21 20:34:58  wasistho
! added new initflo parameter
!
! Revision 1.3  2005/09/09 03:29:44  wasistho
! added INITFLO_PISTON_EXPAN
!
! Revision 1.2  2005/05/02 18:09:34  wasistho
! added cylindrical Taylor inflow profile capability
!
! Revision 1.1  2005/04/29 03:32:33  wasistho
! added distribution bc file generator
!
!
!
!******************************************************************************







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
! Purpose: define various paramters pertinent to PERI
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: PERI_ModParameters.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE PERI_ModParameters

  USE ModDataTypes
  IMPLICIT NONE

! PERI integer parameters -----------------------------------------------------

! general:  

  INTEGER, PARAMETER :: PERI_FLOW_NONE    = 0, &     ! non active
                        PERI_FLOW_CPR     = 1, &     ! CPR flow
                        PERI_FLOW_CHANNEL = 2, &     ! turb. channel flow
                        PERI_FLOW_BOLA    = 3        ! boundary layer flow
! specific to cpr:

  INTEGER, PARAMETER :: CPR_RHO           = 1,  &
                        CPR_RUC           = 2,  & 
                        CPR_RVC           = 3,  & 
                        CPR_UVE           = 4,  & 
                        CPR_VVE           = 5,  & 
                        CPR_TMP           = 6,  & 
                        CPR_PRS           = 7,  & 
                        CPR_DOR           = 8,  & 
                        CPR_DOU           = 9,  & 
                        CPR_DOV           = 10, & 
                        CPR_DOT           = 11, & 
                        CPR_NVAR          = 7,  &    ! CPR number of variables
                        CPR_NCOMP         = 11       ! CPR total variables

  INTEGER, PARAMETER :: GAS_NVAR          = 2        ! # snd/rcv gas variables

! specific to channel:

  INTEGER, PARAMETER :: CNL_PGRAD_MASSFLX = 0,  &    ! type of pgrad calc.meth.
                        CNL_PGRAD_TAUWALL = 1   

! PERI real parameters --------------------------------------------------------

! general:
  REAL(RFREAL), PARAMETER :: PERI_REAL_SMALL = 1.E-16_RFREAL ! small real nmbr

! specific to channel flow:
  REAL(RFREAL), PARAMETER :: CNL_CRITREYN = 1150._RFREAL     ! critical Re

! file IDs --------------------------------------------------------------------

  INTEGER, PARAMETER :: IF_PERI           = 100

END MODULE PERI_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_ModParameters.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2003/09/18 01:56:38  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************







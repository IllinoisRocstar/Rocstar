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
! Purpose: Initialize user input and parameters derived from it for RADI to 
!          default values.
!
! Description: User input and derived parameters are set to default before 
!              overruled by user input
!
! Input: regions data
!
! Output: regions = initial/default values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_InitInputValues.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_InitInputValues( regions )
#endif
#ifdef RFLU
SUBROUTINE RADI_InitInputValues
#endif

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
  USE ModRadiation, ONLY  : t_radi_input
#endif
#ifdef RFLU
  USE ModRadiation, ONLY  : radiInput
#endif
  USE ModError
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region), POINTER :: regions(:)
#endif

! ... loop variables
#ifdef RFLO
  INTEGER :: iReg
#endif
  INTEGER :: m, n

! ... local variables
  TYPE(t_global), POINTER     :: global
  TYPE(t_radi_input), POINTER :: input

  INTEGER :: errorFlag

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RADI_InitInputValues',&
  'RADI_InitInputValues.F90' )

! global values ---------------------------------------------------------------

  global%radiActive = .FALSE.

! region related values -------------------------------------------------------

#ifdef RFLO
  DO iReg=1,global%nRegions

    input => regions(iReg)%radiInput
#endif
#ifdef RFLU
    input => radInput
#endif

! - general

    input%radiModel  = RADI_MODEL_NONE
    input%media      = RADI_MEDIA_ARTIF

! - flux limited diffusion parameter values

    input%fluxLim    = FLD_LIM_LP
    input%spaceDiscr = FLD_DISCR_CEN
    input%spaceOrder = FLD_DISCR_ORD1
    input%vis2       = 0.0_RFREAL
    input%vis4       = 0.0_RFREAL
    input%smoocf     = -1._RFREAL

! - RTE related initial values

    input%solMethod  = RADI_NUM_NONE
    input%nOrdin     = 1
    input%nPol       = 1
    input%nAzi       = 1
    input%nAng       = input%nOrdin

! - optical constants

    ALLOCATE( input%optConst(NPROPERTY,NPHASE), stat = errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    DO n = 1,NPHASE
      DO m = 1,NPROPERTY
        input%optConst(m,n) = 0._RFREAL
        IF (m == PHASE_PROP_D) input%optConst(m,n) = RADI_REAL_MICRON
      ENDDO ! m
    ENDDO   ! n

! - set initial gas volume fraction to 1.
    input%optConst(PHASE_PROP_V,RADI_PHASE_GAS) = 1._RFREAL 

! - define gas ext. efficiency s.t. extinction coeff.= 1.e-4 (repr.for air) 
    input%optConst(PHASE_PROP_Q,RADI_PHASE_GAS) = 1.E-4_RFREAL/1.5_RFREAL* &
                              input%optConst(PHASE_PROP_D,RADI_PHASE_GAS)/ &
                              input%optConst(PHASE_PROP_V,RADI_PHASE_GAS)
#ifdef RFLO
  ENDDO  ! iReg
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_InitInputValues.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:10:30  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.8  2004/09/22 01:30:43  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.7  2004/09/18 17:40:55  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.6  2003/08/01 22:16:20  wasistho
! prepared rocrad for Genx
!
! Revision 1.5  2003/07/30 22:23:18  wasistho
! enter part and smoke data into radiation
!
! Revision 1.4  2003/07/23 03:13:25  wasistho
! cured baby illness
!
! Revision 1.3  2003/07/22 03:03:47  wasistho
! include logical write-parameter
!
! Revision 1.2  2003/07/18 01:38:54  wasistho
! removed bcModel from input data structure
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************








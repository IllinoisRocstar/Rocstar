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
! Purpose: calls routines for initialisation of RADI solutions
!
! Description: Initiate solution variables (qr, radInd, radCoef) to zero.
!              No initial work required for diffusion approximation method of 
!              radiation. For general (RTE) methods, angular discretization
!              is done once in initialisation, while directional weights are 
!              computed for each mesh configuration (once in initialisation and 
!              called again when the grid move).
!
! Input: iReg   = index of current region
!        region = data of current region
!
! Output: for diffusion approximation: none.
!         for general method (RTE): angular mesh and directional weights.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_InitSolution.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_InitSolution( iReg,region )
#endif
#ifdef RFLU
SUBROUTINE RADI_InitSolution
#endif

  USE ModDataTypes
#ifdef RFLO
  USE ModDataStruct, ONLY      : t_region

#include "Indexing.h"
#endif
  USE ModGlobal, ONLY          : t_global
!  USE RADI_ModInterfaces, ONLY : RADI_AngularMesh, RADI_CalcDirWeights
  USE ModError
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: iReg

! ... local variables
  TYPE(t_global), POINTER     :: global

  INTEGER :: iLev, radiModel

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_InitSolution',&
  'RADI_InitSolution.F90' )

! initiate solution, perform initial works and works needed to be done 
! everytime the grid change

#ifdef RFLO
  iLev      = region%currLevel
  radiModel = region%radiInput%radiModel

  region%levels(iLev)%radi%qri     = 0._RFREAL
  region%levels(iLev)%radi%qrj     = 0._RFREAL
  region%levels(iLev)%radi%qrk     = 0._RFREAL
  region%levels(iLev)%radi%radInt  = 0._RFREAL
  region%levels(iLev)%radi%radCoef = 0._RFREAL
  region%levels(iLev)%radi%wvInt   = 0._RFREAL
  region%levels(iLev)%radi%goFact  = 0._RFREAL

  IF ((radiModel == RADI_MODEL_RTEGRAY)  .OR. &
      (radiModel == RADI_MODEL_RTEBAND)) THEN
!    CALL RADI_angularMesh( region )
!    CALL RADI_calcDirWeights( region )
  ENDIF

#endif
#ifdef RFLU
  radiModel = radiInput%radiModel

#endif

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_InitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_InitSolution.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.3  2003/07/30 22:22:53  wasistho
! enter part and smoke data into radiation
!
! Revision 1.2  2003/07/23 03:14:10  wasistho
! cured baby illness
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************








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
! Purpose: Set values derived from user input.
!
! Description: Derived variables are set based on user input parameters.
!              Derived variables are components of radiInput data type 
!              (region%radiInput).
!
! Input: regions = input parameters for all regions.
!
! Output: regions = derived variables stored as part of radiInput data.
!
! Notes: Unlike mixture, derived parameters/variables are not necessarily dv.
!
!******************************************************************************
!
! $Id: RADI_DerivedInputValues.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_DerivedInputValues

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModRadiation, ONLY  : t_radi_input
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... local variables
  TYPE(t_global), POINTER     :: global
  TYPE(t_radi_input), POINTER :: input

  INTEGER      :: errorFlag, angles(3,2)
  REAL(RFREAL) :: pi, twopi

!******************************************************************************

  CALL RegisterFunction( global,'RADI_DerivedInputValues',&
  'RADI_DerivedInputValues.F90' )

! set local constants ---------------------------------------------------------

  pi    = global%pi
  twopi = 2._RFREAL*pi

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__, &
                  'Error in reading real numbers from string' )
20   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__, &
                  'Number of intensity angles is inconsistent' )

999  CONTINUE

END SUBROUTINE RADI_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_DerivedInputValues.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:10:04  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.6  2004/09/22 01:31:23  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.5  2004/09/18 17:41:18  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.4  2003/08/01 22:16:10  wasistho
! prepared rocrad for Genx
!
! Revision 1.3  2003/07/23 03:13:49  wasistho
! cured baby illness
!
! Revision 1.2  2003/07/22 03:05:41  wasistho
! include logical write-parameter
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************








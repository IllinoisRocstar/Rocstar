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
! ******************************************************************************
!
! Purpose: Collection of initialization routines.
!
! Description: global, region and patch data are initialized.
!
! Notes: initialization of conservative variables are done through prep tools.
!
! Todo: to be moved libflo/RFLO_InitInputValues to this module
!
! ******************************************************************************
!
! $Id: RFLO_ModInitialization.F90,v 1.4 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModInitialization

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_InitNonCvSolution

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModInitialization.F90,v $ $Revision: 1.4 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

  
!******************************************************************************
!
! Purpose: Initialisation of non CV solution
!
! Description: none.
!
! Input: region = data of current region and its patches
!
! Output: region, patches = selected region and patch data initialized
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_InitNonCvSolution( region )

  USE RFLO_ModPatchAeroCoeffs, ONLY : RFLO_InitPatchAeroCoeffs

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  INTEGER :: iLev

  TYPE(t_global), POINTER :: global
  TYPE(t_patch),  POINTER :: patch

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLO_ModInitialization.F90,v $ $Revision: 1.4 $'

  global => region%global
  CALL RegisterFunction( global,'RFLO_InitNonCvSolution',&
       'RFLO_ModInitialization.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

! global data -----------------------------------------------------------------

! region data -----------------------------------------------------------------

! patch data ------------------------------------------------------------------

  DO iPatch=1,region%nPatches
    patch => region%levels(iLev)%patches(iPatch)
    CALL RFLO_InitPatchAeroCoeffs( region,patch )
  ENDDO ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_InitNonCvSolution


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModInitialization

! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLO_ModInitialization.F90,v $
!   Revision 1.4  2008/12/06 08:44:16  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:27  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2006/03/24 04:58:45  wasistho
!   modified
!
!   Revision 1.1  2006/03/22 03:06:50  wasistho
!   initial import
!
!
! ******************************************************************************








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
! Purpose: compute RHS for Lagrangian particle position field.
!
! Description: none.
!
! Input: region  = data of current region.
!
! Output: region%levels%plag%rhs
!
! Notes: Use negative values of rhs for consistent RK Updating step.
!
!******************************************************************************
!
! $Id: PLAG_CalcRhsPosition.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_calcRhsPosition( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER :: nPcls
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pDv, pRhs
  
  TYPE(t_plag) ,  POINTER :: pPlag  
  TYPE(t_global), POINTER :: global 
   
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcRhsPosition.F90,v $ $Revision: 1.3 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_calcRhsPosition',&
  'PLAG_CalcRhsPosition.F90' )

! Get dimensions --------------------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel    
  nPcls = region%levels(iLev)%plag%nPcls
#endif
#ifdef RFLU
  nPcls = region%plag%nPcls
#endif

! Set pointers ----------------------------------------------------------------

#ifdef RFLO    
  pPlag => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pPlag => region%plag
#endif
  pDv   => pPlag%dv 
  pRhs  => pPlag%rhs

! Calculate RHS for position vector field -------------------------------------

  DO iPcls = 1, nPcls      
    pRhs(CV_PLAG_XPOS,iPcls) = -pDv(DV_PLAG_UVEL,iPcls)    
    pRhs(CV_PLAG_YPOS,iPcls) = -pDv(DV_PLAG_VVEL,iPcls)       
    pRhs(CV_PLAG_ZPOS,iPcls) = -pDv(DV_PLAG_WVEL,iPcls)         
  END DO  ! iPcls

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_calcRhsPosition

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcRhsPosition.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:03  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.1  2002/10/25 14:14:44  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








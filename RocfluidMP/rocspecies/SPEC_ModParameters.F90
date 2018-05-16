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
! Purpose: Define parameters for species.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_ModParameters.F90,v 1.3 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE SPEC_ModParameters

  IMPLICIT NONE
 
  INTEGER, PARAMETER :: EEV_SPEC_XVEL = 1, &        
                        EEV_SPEC_YVEL = 2, &               
                        EEV_SPEC_ZVEL = 3, &               
                        EEV_SPEC_TEMP = 4, & 
                        EEV_SPEC_NVAR = 4

END MODULE SPEC_ModParameters

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ModParameters.F90,v $
! Revision 1.3  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/11/27 01:47:26  haselbac
! Initial revision
!
! ******************************************************************************







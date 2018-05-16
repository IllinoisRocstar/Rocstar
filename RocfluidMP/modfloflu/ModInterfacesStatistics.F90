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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesStatistics.F90,v 1.6 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesStatistics

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! time statistics
! =============================================================================

#ifdef STATS
  SUBROUTINE ReadStatisticSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadStatisticSection

#ifdef GENX
  SUBROUTINE GenxStatNaming( global, fluidType )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    INTEGER :: fluidType
  END SUBROUTINE GenxStatNaming
#endif
#endif

  END INTERFACE

END MODULE ModInterfacesStatistics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesStatistics.F90,v $
! Revision 1.6  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/12/28 20:30:18  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.3  2004/12/01 00:08:59  wasistho
! added BuildVersionString
!
! Revision 1.2  2004/06/07 23:07:48  wasistho
! added genxStatNaming interface
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************







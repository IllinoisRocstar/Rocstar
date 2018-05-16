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
! Purpose: Write total mass and related info to file.
!
! Description: None.
!
! Input: 
!   regions     Region data
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: WriteTotalMass.F90,v 1.4 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteTotalMass(regions)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModMPI
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  IMPLICIT NONE

! ... parameters

  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'WriteTotalMass',&
  'WriteTotalMass.F90')

! steady flow -----------------------------------------------------------------

  IF ( global%flowType == FLOW_STEADY .AND. & 
       global%myProcid == MASTERPROC ) THEN
    WRITE(IF_MASS,'(I6,4(1X,E23.16))') global%currentIter,global%totalMass, &
                                       global%massIn,global%massOut, & 
                                       global%totalVol

! unsteady flow ---------------------------------------------------------------

  ELSE IF ( global%flowType == FLOW_UNSTEADY .AND. & 
            global%myProcid == MASTERPROC ) THEN
    WRITE(IF_MASS,'(5(1X,E23.16))') global%currentTime,global%totalMass, &
                                    global%massIn,global%massOut, & 
                                    global%totalVol
  END IF ! global%flowType

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE WriteTotalMass

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteTotalMass.F90,v $
! Revision 1.4  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:52:30  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.3  2002/11/15 21:29:33  haselbac
! Deleted RFLU stuff (moved elsewhere), now only write
!
! Revision 1.2  2002/11/15 14:09:25  haselbac
! Changed output format
!
! Revision 1.1  2002/11/08 21:55:48  haselbac
! Initial revision
!
!******************************************************************************








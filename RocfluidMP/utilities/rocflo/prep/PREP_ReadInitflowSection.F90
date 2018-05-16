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
! Purpose: read in user input related to flow initialization.
!
! Description: none.
!
! Input: regions
!
! Output: regions
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_ReadInitflowSection.F90,v 1.7 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadInitflowSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE PREP_ModInterfaces, ONLY : ReadregionSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(10) :: keys(13)

  INTEGER :: brbeg, brend

  LOGICAL :: defined(13)

  REAL(RFREAL) :: vals(13)

!******************************************************************************

  CALL RegisterFunction( regions(1)%global,'ReadInitflowSection', &
                         'PREP_ReadInitflowSection.F90' )

! specify keywords and search for them

  keys(1)  = 'NDUMMY'
  keys(2)  = 'VELX'
  keys(3)  = 'VELY'
  keys(4)  = 'VELZ'
  keys(5)  = 'PRESS'
  keys(6)  = 'DENS'
  keys(7)  = 'XSPLIT'
  keys(8)  = 'DSVELX'
  keys(9)  = 'DSVELY'
  keys(10) = 'DSVELZ'
  keys(11) = 'DSPRESS'
  keys(12) = 'DSDENS'
  keys(13) = 'FLOWCASE'

  CALL ReadregionSection( regions(1)%global,IF_INPUT,13,keys,vals, &
                          brbeg,brend,defined )

  regions(brbeg:brend)%mixtInput%iniXsplit = HUGE(1._RFREAL)

  IF (defined(1)) regions(brbeg:brend)%nDumCells = ABS(vals(1)+0.5_RFREAL)
  IF (defined(2)) regions(brbeg:brend)%mixtInput%iniVelX   = vals(2)
  IF (defined(3)) regions(brbeg:brend)%mixtInput%iniVelY   = vals(3)
  IF (defined(4)) regions(brbeg:brend)%mixtInput%iniVelZ   = vals(4)
  IF (defined(5)) regions(brbeg:brend)%mixtInput%iniPress  = vals(5)
  IF (defined(6)) regions(brbeg:brend)%mixtInput%iniDens   = vals(6)
  IF (defined(7)) regions(brbeg:brend)%mixtInput%iniXsplit = vals(7)
  IF (defined(8)) regions(brbeg:brend)%mixtInput%iniVelX2  = vals(8)
  IF (defined(9)) regions(brbeg:brend)%mixtInput%iniVelY2  = vals(9)
  IF (defined(10)) regions(brbeg:brend)%mixtInput%iniVelZ2 = vals(10)
  IF (defined(11)) regions(brbeg:brend)%mixtInput%iniPress2= vals(11)
  IF (defined(12)) regions(brbeg:brend)%mixtInput%iniDens2 = vals(12)
  IF (defined(13)) regions(brbeg:brend)%mixtInput%prepIniCase  = &
                                                 ABS(vals(13)+0.5_RFREAL)

! finalize

  CALL DeregisterFunction( regions(1)%global )

END SUBROUTINE ReadInitFlowSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_ReadInitflowSection.F90,v $
! Revision 1.7  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/09/09 03:28:52  wasistho
! added FLOWCASE
!
! Revision 1.4  2005/08/17 21:57:28  wasistho
! changed velx2,etc to dsvelx,etc for safety
!
! Revision 1.3  2005/08/03 17:51:48  wasistho
! added domain-splitting feature to initial condition
!
! Revision 1.2  2004/12/03 03:30:14  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.7  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.6  2003/03/20 22:27:56  haselbac
! Renamed ModInterfaces
!
! Revision 1.5  2003/03/20 19:44:22  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.4  2003/03/20 19:35:43  haselbac
! Modified RegFun call to avoid probs with long 'PREP_ReadInitflowSection.F90' names
!
! Revision 1.3  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/02 15:57:08  jblazek
! Added flow initialization.
!
!******************************************************************************









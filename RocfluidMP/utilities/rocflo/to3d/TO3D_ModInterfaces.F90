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
! $Id: TO3D_ModInterfaces.F90,v 1.4 2008/12/06 08:44:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE TO3D_ModInterfaces

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString

  SUBROUTINE RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                      kpnbeg,kpnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensPhysNodes

  SUBROUTINE RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iLev, iNodeOffset, ijNodeOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetNodeOffset

  SUBROUTINE RFLO_WriteGridRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteGridRegion

  END INTERFACE

END MODULE TO3D_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TO3D_ModInterfaces.F90,v $
! Revision 1.4  2008/12/06 08:44:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:37:59  wasistho
! rflo_modinterfacesto3d to to3d_modinterfaces
!
! Revision 1.1  2004/12/03 02:50:44  wasistho
! added prefix
!
! Revision 1.1  2003/03/20 22:36:54  haselbac
! Initial revision
!
! Revision 1.4  2002/06/14 17:46:51  jblazek
! Added version string.
!
! Revision 1.3  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.2  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2001/12/21 23:56:52  jblazek
! Added utility to convert 2D grids to 3D.
!
!******************************************************************************







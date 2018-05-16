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
! $Id: POST_ModInterfaces.F90,v 1.4 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE POST_ModInterfaces

  IMPLICIT NONE

  INTERFACE

  DOUBLE PRECISION FUNCTION Aver1D( cell,var )
    USE ModDataTypes
    INTEGER               :: cell(8)
    REAL(RFREAL), POINTER :: var(:)
  END FUNCTION Aver1D

  DOUBLE PRECISION FUNCTION Aver( cell,iEq,var )
    USE ModDataTypes
    INTEGER               :: cell(8), iEq
    REAL(RFREAL), POINTER :: var(:,:)
  END FUNCTION Aver

  DOUBLE PRECISION FUNCTION AverDiv( cell,iEq1,var1,iEq2,var2 )
    USE ModDataTypes
    INTEGER               :: cell(8), iEq1, iEq2
    REAL(RFREAL), POINTER :: var1(:,:), var2(:,:)
  END FUNCTION AverDiv

  SUBROUTINE BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString

  SUBROUTINE RFLO_CopyGeometryDummy( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CopyGeometryDummy

  SUBROUTINE RFLO_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_DerivedInputValues

  SUBROUTINE RFLO_GenerateCoarseGrids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GenerateCoarseGrids

  SUBROUTINE RFLO_GetCellOffset( region,iLev,iCellOffset,ijCellOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iCellOffset, ijCellOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetCellOffset

  SUBROUTINE RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iNodeOffset, ijNodeOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetNodeOffset

  SUBROUTINE RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                                  kdcbeg,kdcend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummy

  SUBROUTINE RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend,jdnbeg,jdnend, &
                                       kdnbeg,kdnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummyNodes

  SUBROUTINE RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                      kpnbeg,kpnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensPhysNodes

  SUBROUTINE RFLO_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_InitInputValues

  SUBROUTINE MixtureProperties( region,inBeg,inEnd,gasUpdate )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: inBeg, inEnd
    LOGICAL :: gasUpdate
    TYPE(t_region) :: region
  END SUBROUTINE MixtureProperties

  SUBROUTINE RFLO_ReadRegionTopology( global,regions )
    USE ModDataStruct, ONLY : t_region
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadRegionTopology

  SUBROUTINE RFLO_ReadGridRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadGridRegion

  SUBROUTINE PrintPostInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PrintPostInput

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE RFLO_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadSolutionRegion

  SUBROUTINE RFLO_ReadStatRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadStatRegion

  SUBROUTINE WriteGeneric( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE WriteGeneric

  SUBROUTINE WriteTecplotAscii( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE WriteTecplotAscii

  SUBROUTINE WriteTecplotBinary( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE WriteTecplotBinary

#ifdef TURB
  SUBROUTINE TURB_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_ReadInputFile

  SUBROUTINE TURB_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_DerivedInputValues

  SUBROUTINE TURB_RansSAGetEddyVis( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansSAGetEddyVis

  SUBROUTINE TURB_RFLO_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_RFLO_ReadSolutionRegion
#endif
  END INTERFACE

END MODULE POST_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_ModInterfaces.F90,v $
! Revision 1.4  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:20:14  wasistho
! rflo_modinterfacespost to post_modinterfaces
!
! Revision 1.1  2004/12/03 02:03:16  wasistho
! added prefix
!
! Revision 1.8  2004/11/09 10:49:49  wasistho
! added statistics to rflopost
!
! Revision 1.7  2004/07/28 01:50:33  wasistho
! added print input
!
! Revision 1.6  2004/07/24 03:50:23  wasistho
! use postSection instead of command line input
!
! Revision 1.5  2004/03/11 03:36:45  wasistho
! changed rocturb nomenclature
!
! Revision 1.4  2004/02/11 03:25:53  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.3  2004/02/07 01:18:29  wasistho
! added turbulence related results in rocflo post processing
!
! Revision 1.2  2003/09/19 22:38:11  jblazek
! Added turbulence input.
!
! Revision 1.1  2003/03/20 22:25:02  haselbac
! Initial revision
!
! Revision 1.11  2002/10/19 00:40:31  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.10  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.8  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.7  2002/07/20 00:42:05  jblazek
! Added ASCII Tecplot format.
!
! Revision 1.6  2002/06/14 17:16:41  jblazek
! Added version string.
!
! Revision 1.5  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.4  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.2  2002/01/12 00:02:49  jblazek
! Added postprocessor.
!
!******************************************************************************







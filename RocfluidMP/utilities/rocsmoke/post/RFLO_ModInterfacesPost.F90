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
! $Id: RFLO_ModInterfacesPost.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************
  
MODULE RFLO_ModInterfacesPost

  IMPLICIT NONE

  INTERFACE

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

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE RFLO_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadSolutionRegion

  SUBROUTINE WriteTecplotAscii( iReg,iLev,plotType,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg, iLev, plotType
    TYPE(t_region) :: region
  END SUBROUTINE WriteTecplotAscii

  SUBROUTINE PEUL_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER,        INTENT(IN) :: iReg
    TYPE(t_region), POINTER    :: regions(:)
  END SUBROUTINE PEUL_ReadSolutionRegion

#ifdef TURB
  SUBROUTINE TURB_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_ReadInputFile

  SUBROUTINE TURB_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_DerivedInputValues
#endif

  END INTERFACE

END MODULE RFLO_ModInterfacesPost

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModInterfacesPost.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2003/09/25 15:40:22  jferry
! Implented Rocsmoke post-processing
!
!
!******************************************************************************







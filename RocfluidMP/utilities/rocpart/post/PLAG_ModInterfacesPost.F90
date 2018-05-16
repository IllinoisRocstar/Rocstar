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
! $Id: PLAG_ModInterfacesPost.F90,v 1.8 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************
  
MODULE PLAG_ModInterfacesPost

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

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE PLAG_BinSortNozzleInlet( iReg,iLev,region,iRegBin )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg, iLev, iRegBin
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_BinSortNozzleInlet

  SUBROUTINE PLAG_BinSortSpatialDist( iReg,iLev,region,iRegBin )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg, iLev, iRegBin
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_BinSortSpatialDist

  SUBROUTINE PLAG_AllocateMemoryPost( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER, INTENT(IN) :: iReg
  END SUBROUTINE PLAG_AllocateMemoryPost
  
  SUBROUTINE PLAG_CalcDerivedVariables( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_CalcDerivedVariables

  SUBROUTINE PLAG_DeallocateMemoryPost( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER, INTENT(IN) :: iReg
  END SUBROUTINE PLAG_DeallocateMemoryPost
  
  SUBROUTINE PLAG_IntrpMixtProperties( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_IntrpMixtProperties

  SUBROUTINE PLAG_ProcessEulerField( regions,iReg,nPclsSum )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN) :: iReg,nPclsSum
  END SUBROUTINE PLAG_ProcessEulerField
  
  SUBROUTINE PLAG_ReadSolutionFilePost( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_ReadSolutionFilePost

#ifdef STATS
  SUBROUTINE PLAG_ReadStatPost( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN) :: iReg
  END SUBROUTINE PLAG_ReadStatPost
#endif

  SUBROUTINE PLAG_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_UserInput

#ifdef STATS
  SUBROUTINE PLAG_WriteStatTecAscii( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_WriteStatTecAscii
#endif
  
  SUBROUTINE PLAG_WriteTecplotAscii( iReg,iLev,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg, iLev
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_WriteTecplotAscii

! *****************************************************************************
!  RFLO-specific routines
! *****************************************************************************

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

  SUBROUTINE RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                                 kpcbeg,kpcend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensPhys

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
  END INTERFACE

END MODULE PLAG_ModInterfacesPost

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModInterfacesPost.F90,v $
! Revision 1.8  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2005/02/16 23:44:16  fnajjar
! Moved statistics-specific routines inside ifdef construct
!
! Revision 1.5  2005/02/16 14:51:07  fnajjar
! Added interface calls to Tecplot-based statistics file
!
! Revision 1.4  2004/11/17 22:14:45  fnajjar
! Added interface calls to eulerian capability
!
! Revision 1.3  2004/11/13 21:59:35  fnajjar
! Added interface for for PLAG_BinSortSpatialDist
!
! Revision 1.2  2004/05/24 14:24:53  fnajjar
! Included interface and call to binning routine
!
! Revision 1.1.1.1  2003/05/06 16:14:38  fnajjar
! Import of postprocessing tool for Rocpart
!
!******************************************************************************







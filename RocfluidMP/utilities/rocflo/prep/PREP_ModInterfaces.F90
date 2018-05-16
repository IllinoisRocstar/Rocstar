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
! $Id: PREP_ModInterfaces.F90,v 1.6 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE PREP_ModInterfaces

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString

  SUBROUTINE CheckBcValidity( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE CheckBcValidity

  SUBROUTINE GetGrid( gridLevel,iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: gridLevel,iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE GetGrid

  SUBROUTINE InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE InitInputValues

  SUBROUTINE InitializeFlowField( iLev,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iLev
    TYPE(t_region) :: region
  END SUBROUTINE InitializeFlowField

  SUBROUTINE PrintPrepInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PrintPrepInput

  SUBROUTINE ReadBcInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadBcInputFile

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE ReadFormatsSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadFormatsSection

  SUBROUTINE ReadInitflowSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInitflowSection

  SUBROUTINE ReadMultigridSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadMultigridSection

  SUBROUTINE ReadReferenceSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadReferenceSection

  SUBROUTINE ReadTimestepSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTimestepSection

  SUBROUTINE ReadRegionSection( global,fileID,nvals,keys,vals, &
                                brbeg,brend,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, brbeg, brend
    CHARACTER(*) :: keys(nvals)
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadRegionSection

  SUBROUTINE ReadSection( global,fileID,nvals,keys,vals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals
    CHARACTER(*) :: keys(nvals)
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadSection

  SUBROUTINE RFLO_CopyGeometryDummy( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CopyGeometryDummy

  SUBROUTINE RFLO_ExtrapolateGeometry( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_ExtrapolateGeometry

  SUBROUTINE RFLO_GenerateCoarseGrids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GenerateCoarseGrids

  SUBROUTINE RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                                  kdcbeg,kdcend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummy

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

  SUBROUTINE RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend,jdnbeg,jdnend,&
                                       kdnbeg,kdnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummyNodes

  SUBROUTINE RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                        jbeg,jend,kbeg,kend )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_GetPatchIndicesNodes

  SUBROUTINE RFLO_ReadGridRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadGridRegion

  SUBROUTINE RFLO_ReadRegionTopology( global,regions )
    USE ModDataStruct, ONLY : t_region
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadRegionTopology

  SUBROUTINE RFLO_WriteSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteSolutionRegion

#ifdef GENX
  SUBROUTINE RFLO_InitGenxInterfacePrep( iReg,region,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region) :: region
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE RFLO_InitGenxInterfacePrep 
                                   
  SUBROUTINE GenxInitSolution( gridLevel,iReg,regions,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    INTEGER :: gridLevel, iReg
    TYPE(t_region), POINTER :: regions(:)
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE GenxInitSolution

  SUBROUTINE GenxWriteSolution( gridLevel,iReg,region,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    INTEGER :: gridLevel, iReg
    TYPE(t_region) :: region
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE GenxWriteSolution

  SUBROUTINE GenxWriteRocinout( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE GenxWriteRocinout
#endif

  END INTERFACE

END MODULE PREP_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_ModInterfaces.F90,v $
! Revision 1.6  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/09/30 00:14:24  wasistho
! added post read grid processing routines
!
! Revision 1.3  2005/05/02 18:07:45  wasistho
! added cylindrical Taylor inflow profile capability
!
! Revision 1.2  2004/12/03 03:29:57  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.5  2004/07/27 20:28:48  wasistho
! added readBcInputFile and checkBcValidity
!
! Revision 1.4  2004/07/27 03:36:35  wasistho
! add printPrepInput
!
! Revision 1.3  2004/07/23 04:31:18  wasistho
! Genx: readin from Rocin, standalone: read .inp file i.o. command line input
!
! Revision 1.2  2004/06/30 00:00:12  wasistho
! migrated to Roccom-3
!
! Revision 1.1  2003/03/20 22:29:29  haselbac
! Initial revision
!
! Revision 1.5  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/06/14 17:35:11  jblazek
! Added version string.
!
! Revision 1.3  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/02 15:57:08  jblazek
! Added flow initialization.
!
!******************************************************************************







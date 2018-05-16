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
! $Id: TFLU_ModInterfaces.F90,v 1.5 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE TFLU_ModInterfaces

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString

  SUBROUTINE ConvertFlo2FluMesh( iFlag,iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iFlag, iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ConvertFlo2FluMesh

  SUBROUTINE ConvertFlo2FluPatch( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ConvertFlo2FluPatch

  SUBROUTINE CorrectNedges( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE CorrectNedges

  SUBROUTINE GetBndVertType( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE GetBndVertType

  SUBROUTINE PrintTofluInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PrintTofluInput

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE ReadFormatsSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadFormatsSection

  SUBROUTINE ReadMultigridSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadMultigridSection

  SUBROUTINE RFLO_ReadRegionTopology( global,regions )
    USE ModDataStruct, ONLY : t_region
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadRegionTopology

  SUBROUTINE RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                      kpnbeg,kpnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensPhysNodes

  SUBROUTINE RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend,jdnbeg,jdnend,&
                                       kdnbeg,kdnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummyNodes

  SUBROUTINE RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iNodeOffset, ijNodeOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetNodeOffset

  SUBROUTINE RFLO_GetPatchDirection( patch,idir,jdir,kdir )
    USE ModBndPatch, ONLY : t_patch
    INTEGER       :: idir, jdir, kdir
    TYPE(t_patch) :: patch
  END SUBROUTINE RFLO_GetPatchDirection

  SUBROUTINE RFLO_GetPatchIndices( region,patch,iLev, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_GetPatchIndices

  SUBROUTINE RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                        ibeg,iend,jbeg,jend,kbeg,kend )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_GetPatchIndicesNodes

  SUBROUTINE RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                                   idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                                   ibeg,iend,jbeg,jend,kbeg,kend, &
                                   ibegSrc,iendSrc,jbegSrc,jendSrc, &
                                   kbegSrc,kendSrc,mapMat )
    INTEGER :: lb, lbs, l1SrcDir, l2SrcDir
    INTEGER :: idir, jdir, kdir, idirSrc, jdirSrc, kdirSrc
    INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
    INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc
    INTEGER :: mapMat(3,4)
    LOGICAL :: align
  END SUBROUTINE RFLO_GetPatchMapping

  SUBROUTINE RFLO_ReadGridRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadGridRegion

  SUBROUTINE WriteFluCellMap( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteFluCellMap

  SUBROUTINE WriteFluDimens( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteFluDimens

  SUBROUTINE WriteFluGrid( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteFluGrid

  END INTERFACE

END MODULE TFLU_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_ModInterfaces.F90,v $
! Revision 1.5  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/12/21 22:38:11  wasistho
! added writeFluCellMap
!
! Revision 1.2  2004/12/03 03:44:15  wasistho
! rflo_modinterfacestoflu to tflu_modinterfaces
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.2  2004/08/18 02:10:09  wasistho
! added new routines to create dimension file
!
! Revision 1.1.1.1  2004/08/17 01:41:39  wasistho
! initial checkin
!
!
!******************************************************************************







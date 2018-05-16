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
! $Id: SURF_ModInterfaces.F90,v 1.4 2008/12/06 08:44:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE SURF_ModInterfaces

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString

  SUBROUTINE RFLO_CopyGeometryDummy( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CopyGeometryDummy

  SUBROUTINE RFLO_CopyTopologyLevels( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CopyTopologyLevels

  SUBROUTINE CountInteractingPatches( regions,nInteract )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: nInteract
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE CountInteractingPatches

  SUBROUTINE RFLO_GenerateCoarseGrids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GenerateCoarseGrids

  SUBROUTINE RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iNodeOffset, ijNodeOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetNodeOffset

  SUBROUTINE RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend,jdnbeg,jdnend, &
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

  SUBROUTINE RFLO_ReadBcInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcInputFile

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

  SUBROUTINE WriteSurfaceGrid( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE WriteSurfaceGrid

  END INTERFACE

END MODULE SURF_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SURF_ModInterfaces.F90,v $
! Revision 1.4  2008/12/06 08:44:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:35:48  wasistho
! rflo_modinterfacessurf to surf_modinterfaces
!
! Revision 1.1  2004/12/03 02:47:00  wasistho
! added prefix
!
! Revision 1.1  2003/03/20 22:34:11  haselbac
! Initial revision
!
! Revision 1.1  2002/10/19 00:40:31  jblazek
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







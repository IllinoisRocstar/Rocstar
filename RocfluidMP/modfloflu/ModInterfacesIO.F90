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
! Purpose: set explicit interfaces to subroutines and functions
!          related to I/O.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesIO.F90,v 1.9 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesIO

  IMPLICIT NONE

  INTERFACE

  INTEGER FUNCTION BuildPatchIdentifier(iRegion,iPatch)
    INTEGER, INTENT(IN) :: iPatch,iRegion
  END FUNCTION BuildPatchIdentifier

  SUBROUTINE MakeNumberedKeys(keys,indBegin,string,numBegin,numEnd,numSkip)
    CHARACTER(*), INTENT(inout) :: keys(:)
    CHARACTER(*), INTENT(in)    :: string
    INTEGER,      INTENT(in)    :: indBegin,numBegin,numEnd,numSkip
  END SUBROUTINE MakeNumberedKeys

  SUBROUTINE ReadBothRegionSection( global,fileID,nvals,nStrVals,keys, &
                                    strKeys,vals,strVals,brbeg,brend,  &
                                    defined,strDefined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, nStrVals, brbeg, brend
    CHARACTER(*) :: keys(nvals), strKeys(nStrVals)
    LOGICAL      :: defined(nvals), strDefined(nStrVals)
    REAL(RFREAL) :: vals(nvals)
    CHARACTER(*) :: strVals(nStrVals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadBothRegionSection

  SUBROUTINE ReadBothSection( global,fileID,nvals,nStrVals,keys,strKeys, &
                              vals,strVals,defined,strDefined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, nStrVals
    CHARACTER(*) :: keys(nvals), strKeys(nStrVals)
    LOGICAL      :: defined(nvals), strDefined(nStrVals)
    REAL(RFREAL) :: vals(nvals)
    CHARACTER(*) :: strVals(nStrVals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadBothSection

  SUBROUTINE ReadAccelerationSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadAccelerationSection

  SUBROUTINE ReadFlowSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadFlowSection

  SUBROUTINE ReadForcesSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadForcesSection

  SUBROUTINE ReadFormatsSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadFormatsSection

#ifdef RFLO
  SUBROUTINE ReadGridMotionSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadGridMotionSection
#endif
#ifdef RFLU
  SUBROUTINE ReadGridMotionSection(regions)
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadGridMotionSection
#endif

  SUBROUTINE ReadInitFlowSection(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions
  END SUBROUTINE ReadInitFlowSection

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE ReadListSection( global,fileID,key,nCols,nRows,vals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nCols, nRows
    CHARACTER(*) :: key
    LOGICAL      :: defined
    REAL(RFREAL), POINTER :: vals(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadListSection
  
  SUBROUTINE ReadMixtureSection(regions)
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadMixtureSection

  SUBROUTINE ReadMultigridSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadMultigridSection

  SUBROUTINE ReadNumericsSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadNumericsSection

#ifdef RFLO
  SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals,brbeg,brend, &
                               prbeg,prend,distrib,fname,defined )
#endif
#ifdef RFLU
  SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals, &
                               prbeg,prend,distrib,fname,bcName,defined )
#endif
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
#ifdef RFLO
    INTEGER      :: brbeg, brend
#endif    
    INTEGER      :: fileID, nvals, prbeg, prend, distrib
    CHARACTER(*) :: keys(nvals), fname
#ifdef RFLU
    CHARACTER(*) :: bcName
#endif    
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPatchSection

  SUBROUTINE ReadPostSection(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPostSection

  SUBROUTINE ReadPrefixedListSection( global,fileID,key,nCols,nRows, &
                                      vals,strVals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nCols, nRows
    CHARACTER(*) :: key
    LOGICAL      :: defined
    REAL(RFREAL), POINTER :: vals(:,:)
    CHARACTER(*), POINTER :: strVals(:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPrefixedListSection

  SUBROUTINE ReadPrepSection(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPrepSection

  SUBROUTINE ReadProbeSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadProbeSection

  SUBROUTINE ReadRandomSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadRandomSection

  SUBROUTINE ReadReferenceSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadReferenceSection

#ifdef RFLU
  SUBROUTINE ReadTimeZoomingSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTimeZoomingSection

  SUBROUTINE ReadRocketSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadRocketSection
#endif

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

  SUBROUTINE ReadStringSection( global,fileID,nvals,keys,vals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER        :: fileID, nvals
    CHARACTER(*)   :: keys(nvals)
    LOGICAL        :: defined(nvals)
    CHARACTER(*)   :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadStringSection

  SUBROUTINE ReadThrustSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadThrustSection

  SUBROUTINE ReadTimestepSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTimestepSection

  SUBROUTINE ReadTransformSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTransformSection

  SUBROUTINE ReadViscositySection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadViscositySection

  SUBROUTINE WriteConvergence( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteConvergence

  SUBROUTINE WriteTotalMass(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE WriteTotalMass

  SUBROUTINE WriteProbe( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE WriteProbe

  SUBROUTINE WriteThrust( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteThrust

  END INTERFACE

END MODULE ModInterfacesIO

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesIO.F90,v $
! Revision 1.9  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2007/04/14 14:29:34  mtcampbe
! Updated for TZ and rocket case constraints
!
! Revision 1.6  2005/10/31 19:26:48  haselbac
! Added interface for ReadMixtureSection
!
! Revision 1.5  2004/06/16 20:00:54  haselbac
! Removed buildFileNameXXX routines
!
! Revision 1.4  2004/04/08 03:17:07  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.3  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.2  2003/08/28 20:05:39  jblazek
! Added acceleration terms.
!
! Revision 1.1  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
!******************************************************************************







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
! $Id: RVAV_ModInterfaces.F90,v 1.10 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************
  
MODULE RVAV_ModInterfaces

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE BuildVersionString

  SUBROUTINE RVAV_ReadFileStream1( regions )
    USE      ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regions(:)
  END SUBROUTINE RVAV_ReadFileStream1
  
  SUBROUTINE RVAV_ComputeAnalyticalSolution( similarityType, regions )
    USE      ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regions(:)
    INTEGER, INTENT(IN)        :: similarityType
  END SUBROUTINE RVAV_ComputeAnalyticalSolution
    
  SUBROUTINE RVAV_ReadFileStream2( regionsS1, regionsS2 )
    USE      ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regionsS1(:)
    TYPE(t_region),      POINTER :: regionsS2(:)
  END SUBROUTINE RVAV_ReadFileStream2
    
  SUBROUTINE RVAV_ReadFileStream2Analyt( global, regionsS2 )
    USE      ModDataStruct, ONLY : t_region
    USE      ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER      :: global
    TYPE(t_region),      POINTER :: regionsS2(:)
  END SUBROUTINE RVAV_ReadFileStream2Analyt
        
  SUBROUTINE RVAV_ReadFileStream2Comput( regionsS1, regionsS2 )
    USE      ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regionsS1(:)
    TYPE(t_region),      POINTER :: regionsS2(:)
  END SUBROUTINE RVAV_ReadFileStream2Comput
        
  SUBROUTINE RVAV_ReadFileStream2Experm( global, regionsS2 )
    USE      ModDataStruct, ONLY : t_region
    USE      ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER      :: global
    TYPE(t_region),      POINTER :: regionsS2(:)
  END SUBROUTINE RVAV_ReadFileStream2Experm  
    
  SUBROUTINE RVAV_ExtractVariables( global, region, &
                                    ibeg, iend, ijump, &
                                    jbeg, jend, jjump, &
                                    kbeg, kend, kjump, &
                                    iOffset, ijOffset, &
                                    variableIndex,     &
                                    fileType,          &
                                    indCp, indMol, ev )
  
    USE ModDataTypes
    USE ModDataStruct, ONLY    : t_region
    USE      ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER      :: global
    TYPE(t_region), INTENT(IN) :: region
    INTEGER, INTENT(IN)        :: ibeg, iend, ijump
    INTEGER, INTENT(IN)        :: jbeg, jend, jjump
    INTEGER, INTENT(IN)        :: kbeg, kend, kjump
    INTEGER, INTENT(IN)        :: iOffset, ijOffset
    INTEGER, INTENT(IN)        :: indCp, indMol
    INTEGER, INTENT(IN)        :: variableIndex
    INTEGER, INTENT(IN)        :: fileType
    REAL(RFREAL), POINTER      :: ev(:,:,:)
  END SUBROUTINE RVAV_ExtractVariables

  SUBROUTINE RVAV_ComputeSimilarField( global, iNodes,jNodes,kNodes, &
                                       similarityType,       &
                                       variableIndex,        &
                                       ev)
    USE      ModGlobal, ONLY : t_global
    USE      ModParameters
    USE      RVAV_ModGlobal
    TYPE(t_global), POINTER    :: global
    INTEGER, INTENT(IN)   :: iNodes,jNodes,kNodes
    INTEGER, INTENT(IN)        :: variableIndex
    INTEGER, INTENT(IN)        :: similarityType
    REAL(RFREAL), POINTER      :: ev(:,:,:)
  END SUBROUTINE RVAV_ComputeSimilarField
  
  SUBROUTINE RVAV_ComputeError( global, iCompare, iNodes, jNodes, kNodes)
    USE      ModGlobal, ONLY : t_global
    USE      ModParameters
    USE      RVAV_ModGlobal
    TYPE(t_global), POINTER    :: global
    INTEGER, INTENT(IN)   :: iCompare
    INTEGER, INTENT(IN)   :: iNodes,jNodes,kNodes
  END SUBROUTINE RVAV_ComputeError
  
  SUBROUTINE RVAV_PlotResults( global, iCompare, iNodes, jNodes, kNodes)
    USE      ModGlobal, ONLY : t_global
    USE      ModParameters
    USE      RVAV_ModGlobal
    TYPE(t_global), POINTER  :: global
    INTEGER, INTENT(IN)   :: iCompare
    INTEGER, INTENT(IN)   :: iNodes,jNodes,kNodes
  END SUBROUTINE RVAV_PlotResults
  
  SUBROUTINE RVAV_BlasiusSolution( fname, regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regions(:)
    CHARACTER(*) :: fname
  END SUBROUTINE RVAV_BlasiusSolution  
  
  SUBROUTINE RVAV_GammBumpSolution( fname, regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regions(:)
    CHARACTER(*) :: fname
  END SUBROUTINE RVAV_GammBumpSolution
    
  SUBROUTINE RVAV_ProudmanCulickSolution( fname, regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region),      POINTER :: regions(:)
    CHARACTER(*) :: fname
  END SUBROUTINE RVAV_ProudmanCulickSolution
  
  SUBROUTINE RVAV_ReadInputFile( global )
    USE      ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER    :: global
  END SUBROUTINE RVAV_ReadInputFile
  
  SUBROUTINE RVAV_ReadSectionStream1( global )
    USE      ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER    :: global
  END SUBROUTINE RVAV_ReadSectionStream1
  
  SUBROUTINE RVAV_ReadSectionStream2( global )
    USE      ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER    :: global
  END SUBROUTINE RVAV_ReadSectionStream2
  
  SUBROUTINE RVAV_ReadComparisonsSection( global )
    USE      ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER    :: global
  END SUBROUTINE RVAV_ReadComparisonsSection
              
  END INTERFACE

END MODULE RVAV_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ModInterfaces.F90,v $
! Revision 1.10  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.7  2002/08/15 19:48:06  jblazek
! Implemented grid deformation capability.
!
! Revision 1.6  2002/07/31 02:34:57  f-najjar
! Split Analytical Solutions into individual routines
!
! Revision 1.5  2002/07/16 22:42:19  f-najjar
! Cleanup of ReadFileStream2
!
! Revision 1.4  2002/06/18 03:18:20  f-najjar
! Included RVAV_computeAnalyticalSolution
!
! Revision 1.3  2002/06/15 17:55:54  f-najjar
! Bug fix of Extract Variables and New Call for computSimilarField
!
! Revision 1.2  2002/06/14 17:00:52  jblazek
! Added version string.
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************







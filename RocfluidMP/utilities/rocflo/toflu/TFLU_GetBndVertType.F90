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
! Purpose: get type of boundary vertices, needed to correct global%tofluNedges
!
! Description: the 'type' has 6 components which have initial value of zero. 
!              Component 2, 4, or 6 switches to 1 if the lbound value of the 
!              patch where the boundary vertex resides matches the component 
!              number. This way, the points to be corrected for nEdges can
!              be distinguished. Only non-connecting patches are processed, 
!              as the procedure concerns only the global boundary vertices.
!              The correction of nEdges itself is performed in subroutine
!              CorrectNedges.
!
! Input: iReg    = region number
!        regions = data for all regions
!
! Output: global%tofluBType
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TFLU_GetBndVertType.F90,v 1.4 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GetBndVertType( iReg,Regions )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModBndPatch, ONLY   : t_patch
  USE TFLU_ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch, i, j, k

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_grid) , POINTER  :: grid
  TYPE(t_patch) , POINTER :: patch

  INTEGER :: iLev, lbound, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'GetBndVertType',&
  'TFLU_GetBndVertType.F90' )

! obtain parameters and pointers, and loop over region patches ----------------

  iLev =  regions(iReg)%currLevel
  grid => regions(iReg)%levels(iLev)%grid

  DO iPatch=1,regions(iReg)%nPatches

    patch   => regions(iReg)%levels(iLev)%patches(iPatch)
    lbound  =  patch%lbound
    bcType  =  patch%bcType

! - only at non-connecting patches --------------------------------------------

    IF (bcType<BC_REGIONCONF .OR. bcType>BC_REGIONCONF+BC_RANGE) THEN

      CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                      ibeg,iend,jbeg,jend,kbeg,kend )

! --- identify the vertices where nEdges needs to be corrected 

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend

            IF (lbound==2) THEN
              global%tofluBType(2,grid%tofluLoc2g(i,j,k)) = 1
            ELSEIF (lbound==4) THEN
              global%tofluBType(4,grid%tofluLoc2g(i,j,k)) = 1
            ELSEIF (lbound==6) THEN
              global%tofluBType(6,grid%tofluLoc2g(i,j,k)) = 1
            ENDIF

          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    
    ENDIF  ! bcType
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GetBndVertType

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_GetBndVertType.F90,v $
! Revision 1.4  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:43:58  wasistho
! rflo_modinterfacestoflu to tflu_modinterfaces
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:58:20  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/18 02:15:12  wasistho
! added new routines to create dimension file
!
!
!******************************************************************************








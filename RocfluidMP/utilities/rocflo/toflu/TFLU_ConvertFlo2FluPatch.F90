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
! Purpose: obtain face (quads) and vertex connectivity at non-connecting 
!          boundary patches in structured grid, equivalent to true boundary
!          patches in unstructured grid.
!
! Description: quad and vertex connectivity at boundary patches are obtined
!              from current region l2g mapping
!
! Input: iReg       = region number
!        regions    = data for all regions
!
! Output: global%tofluQuad2v (for quads) and tofluBLoc2g (for vertices)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TFLU_ConvertFlo2FluPatch.F90,v 1.4 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ConvertFlo2FluPatch( iReg,regions )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModBndPatch, ONLY   : t_patch
  USE TFLU_ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndices, &
                                 RFLO_GetPatchIndicesNodes
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

  INTEGER :: iLev, lbound, bcType, iNOff, ijNOff, ibn, ien, ijkN, ng1, ng2
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, igPatch, dims(2), errorFlag
  INTEGER :: iadd, jadd, kadd, iq

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'ConvertFlo2FluPatch',&
  'TFLU_ConvertFlo2FluPatch.F90' )

! obtain boundary vertices and faces (quads) connectivity ---------------------

  iLev = regions(iReg)%currLevel
  grid => regions(iReg)%levels(iLev)%grid

  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

  DO iPatch=1,regions(iReg)%nPatches

    patch   => regions(iReg)%levels(iLev)%patches(iPatch)
    lbound  =  patch%lbound
    bcType  =  patch%bcType

! - only at non-connecting patches --------------------------------------------

    IF (bcType<BC_REGIONCONF .OR. bcType>BC_REGIONCONF+BC_RANGE) THEN

! --- get dimensions needed to output the unstructured grid

      global%tofluNPatches = global%tofluNPatches + 1
      igPatch = global%tofluNPatches 
      CALL RFLO_GetPatchIndices( regions(iReg),patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )

      dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
      dims(2) = ABS(patch%l2end-patch%l2beg) + 2

      global%tofluNbVerts(igPatch) = global%tofluNbVerts(igPatch) + &
                                     dims(1)*dims(2)

      dims(1) = ABS(patch%l1end-patch%l1beg) + 1    ! face values
      dims(2) = ABS(patch%l2end-patch%l2beg) + 1

      global%tofluNbFaces(igPatch) = global%tofluNbFaces(igPatch) + &
                                     dims(1)*dims(2)
      iadd = 0
      jadd = 0
      kadd = 0
      IF (lbound==2) THEN
        iadd=1
      ELSEIF (lbound==4) THEN
        jadd=1
      ELSEIF (lbound==6) THEN
        kadd=1
      ENDIF

! --- obtain quads connectivity depending on lbound ---------------------------

      DO k=kbeg+kadd,kend+kadd
        DO j=jbeg+jadd,jend+jadd
          DO i=ibeg+iadd,iend+iadd

            IF      (lbound==1 .OR. lbound==2) THEN
              ng1 = j - jbeg + 1
              ng2 = k - kbeg + 1
              iq  = (ng2-1)*(jend-jbeg+1)+ng1
              IF (lbound==1) THEN
                global%tofluQuad2v(1,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k  )
                global%tofluQuad2v(2,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k+1)
                global%tofluQuad2v(3,iq,igPatch) = grid%tofluLoc2g(i  ,j+1,k+1)
                global%tofluQuad2v(4,iq,igPatch) = grid%tofluLoc2g(i  ,j+1,k  )
              ELSE
                global%tofluQuad2v(1,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k  )
                global%tofluQuad2v(2,iq,igPatch) = grid%tofluLoc2g(i  ,j+1,k  )
                global%tofluQuad2v(3,iq,igPatch) = grid%tofluLoc2g(i  ,j+1,k+1)
                global%tofluQuad2v(4,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k+1)
              ENDIF
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              ng1 = k - kbeg + 1
              ng2 = i - ibeg + 1
              iq  = (ng2-1)*(kend-kbeg+1)+ng1
              IF (lbound==3) THEN
                global%tofluQuad2v(1,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k  )
                global%tofluQuad2v(2,iq,igPatch) = grid%tofluLoc2g(i+1,j  ,k  )
                global%tofluQuad2v(3,iq,igPatch) = grid%tofluLoc2g(i+1,j  ,k+1)
                global%tofluQuad2v(4,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k+1)
              ELSE
                global%tofluQuad2v(1,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k  )
                global%tofluQuad2v(2,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k+1)
                global%tofluQuad2v(3,iq,igPatch) = grid%tofluLoc2g(i+1,j  ,k+1)
                global%tofluQuad2v(4,iq,igPatch) = grid%tofluLoc2g(i+1,j  ,k  )
              ENDIF
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              ng1 = i - ibeg + 1
              ng2 = j - jbeg + 1
              iq  = (ng2-1)*(iend-ibeg+1)+ng1
              IF (lbound==5) THEN
                global%tofluQuad2v(1,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k  )
                global%tofluQuad2v(2,iq,igPatch) = grid%tofluLoc2g(i  ,j+1,k  )
                global%tofluQuad2v(3,iq,igPatch) = grid%tofluLoc2g(i+1,j+1,k  )
                global%tofluQuad2v(4,iq,igPatch) = grid%tofluLoc2g(i+1,j  ,k  )
              ELSE
                global%tofluQuad2v(1,iq,igPatch) = grid%tofluLoc2g(i  ,j  ,k  )
                global%tofluQuad2v(2,iq,igPatch) = grid%tofluLoc2g(i+1,j  ,k  )
                global%tofluQuad2v(3,iq,igPatch) = grid%tofluLoc2g(i+1,j+1,k  )
                global%tofluQuad2v(4,iq,igPatch) = grid%tofluLoc2g(i  ,j+1,k  )
              ENDIF  ! lbound
            ENDIF    ! lbound
          ENDDO      ! i
        ENDDO        ! j
      ENDDO          ! k

! --- obtain connectivity at boundary vertices --------------------------------

      CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                      ibeg,iend,jbeg,jend,kbeg,kend )
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend

            IF      (lbound==1 .OR. lbound==2) THEN
              ng1 = j - jbeg + 1
              ng2 = k - kbeg + 1
              iq  = (ng2-1)*(jend-jbeg+1)+ng1
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              ng1 = k - kbeg + 1
              ng2 = i - ibeg + 1
              iq  = (ng2-1)*(kend-kbeg+1)+ng1
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              ng1 = i - ibeg + 1
              ng2 = j - jbeg + 1
              iq  = (ng2-1)*(iend-ibeg+1)+ng1
            ENDIF
            global%tofluBLoc2g(iq,igPatch) = grid%tofluLoc2g(i,j,k)

            IF (lbound==2) THEN
              global%tofluMaxBind = &
                  MAX( global%tofluMaxBind , grid%tofluLoc2g(i,j,k) ) 
            ELSEIF (lbound==4) THEN
              global%tofluMaxBind = &
                  MAX( global%tofluMaxBind , grid%tofluLoc2g(i,j,k) ) 
            ELSEIF (lbound==6) THEN
              global%tofluMaxBind = &
                  MAX( global%tofluMaxBind , grid%tofluLoc2g(i,j,k) ) 
            ENDIF
            
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k

      global%tofluIq(igPatch) = iq
    
    ENDIF  ! bcType
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE ConvertFlo2FluPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_ConvertFlo2FluPatch.F90,v $
! Revision 1.4  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:43:49  wasistho
! rflo_modinterfacestoflu to tflu_modinterfaces
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:58:20  wasistho
! lower to upper case
!
! Revision 1.2  2004/08/18 02:10:40  wasistho
! added new routines to create dimension file
!
! Revision 1.1.1.1  2004/08/17 01:41:39  wasistho
! initial checkin
!
!
!******************************************************************************








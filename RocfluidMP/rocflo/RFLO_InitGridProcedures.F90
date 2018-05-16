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
! Purpose: wrap up of procedures involving original initial undeformed grid.
!
! Description: procedures such as computing averaging weights based on initial
!              undeformed grid are called from here.
!
! Input: regions = grid dimensions
!        input from file.
!
! Output: specified in the called routines.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_InitGridProcedures.F90,v 1.8 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_InitGridProcedures( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE RFLO_ModMoveGridFrame, ONLY : RFLO_MgFrameCornPoints, &
                                    RFLO_MgFrameBroadCast, &
                                    RFLO_MgFrameSrchNeighbors 
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_GetGeometry',&
  'RFLO_InitGridProcedures.F90' )

! broadcast block corners and remember closest six neighbors ------------------

  IF (global%moveGridScheme == MOVEGRID_FRAME .OR. &
      global%moveGridScheme == MOVEGRID_FOMS .OR. &
      global%moveGridScheme == MOVEGRID_ELFRAME) THEN
    CALL RFLO_MgFrameCornPoints( regions )
    CALL RFLO_MgFrameBroadCast( regions,0,1 )
    CALL RFLO_MgFrameSrchNeighbors( regions )
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_InitGridProcedures

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InitGridProcedures.F90,v $
! Revision 1.8  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/03/02 01:27:44  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.5  2006/02/08 03:32:17  wasistho
! added movegrid_epde
!
! Revision 1.4  2005/10/28 07:38:54  wasistho
! included moveGridFoms
!
! Revision 1.3  2005/06/10 19:33:51  wasistho
! added call to RFLO_MgFrameCornPoints
!
! Revision 1.2  2005/06/01 07:21:13  wasistho
! set iselect argument in mgFrameBroadcast from 2 to 0
!
! Revision 1.1  2005/05/28 08:10:10  wasistho
! import RFLO_InitGridProcedures
!
!
!******************************************************************************








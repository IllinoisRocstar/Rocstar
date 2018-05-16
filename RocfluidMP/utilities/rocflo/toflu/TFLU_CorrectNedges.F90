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
! Purpose: correct global number of edges
!
! Description: to this point number of edges is calculated by increamenting
!              global%tofluNedges by 3 times number of vertices over all
!              regions subtracted by overcalculated edges due to connecting
!              patches. There remains overcalculated edges at the patches
!              with lbound = 2, 4, and 6. This has to be corrected.
!
! Input: iReg     = region number
!        regions  = data for all regions
!
! Output: updated global%tofluNedges
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TFLU_CorrectNedges.F90,v 1.3 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE CorrectNedges( global )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY     : t_global
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: iBv, lbound

!******************************************************************************

  CALL RegisterFunction( global,'CorrectNedges',&
  'TFLU_CorrectNedges.F90' )

! start -----------------------------------------------------------------------

  DO iBv = 1,global%tofluMaxBind
    DO lbound = 1,6

      IF (lbound==2 .AND. global%tofluBType(2,iBv)==1) THEN
        global%tofluNEdges = global%tofluNEdges - 1

      ELSEIF (lbound==4 .AND. global%tofluBType(4,iBv)==1) THEN
        global%tofluNEdges = global%tofluNEdges - 1

      ELSEIF (lbound==6 .AND. global%tofluBType(6,iBv)==1) THEN
        global%tofluNEdges = global%tofluNEdges - 1

      ENDIF
    ENDDO  ! lbound
  ENDDO    ! iBv

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE CorrectNedges

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_CorrectNedges.F90,v $
! Revision 1.3  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
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








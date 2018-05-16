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
! Purpose: open file for convergence history.
!
! Description: none.
!
! Input: global = case name, steady/unsteady flow.
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_OpenConverFile.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_OpenConverFile( global )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(CHRLEN+4) :: fname

  INTEGER :: errorFlag

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_OpenConverFile',&
  'RFLO_OpenConverFile.F90' )

! open file

  IF (global%myProcid == MASTERPROC) THEN
    fname = TRIM(global%outDir)//TRIM(global%casename)//'.con'

! - append to existing file (restart) or create new file

    IF ((global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL).OR.&
        (global%flowType==FLOW_STEADY   .AND. global%currentIter>1)) THEN
      OPEN(IF_CONVER,file=fname,form='formatted',status='old', &
                     position='append',iostat=errorFlag)
    ELSE
      OPEN(IF_CONVER,file=fname,form='formatted',status='unknown', &
                     iostat=errorFlag)
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_OpenConverFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_OpenConverFile.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.5  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.4  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/02/25 22:36:53  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************








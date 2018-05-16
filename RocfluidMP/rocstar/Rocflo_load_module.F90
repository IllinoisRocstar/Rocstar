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
! Purpose: open windoe to fluids code and register interface functions.
!
! Description: none.
!
! Input: winName = name of fluids window.
!
! Output: registered functions.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Rocflo_load_module.F90,v 1.6 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE Rocflo_load_module( winName )

  USE ModRocstar, ONLY       : t_globalGenx, associate_pointer
  USE ModInterfaces, ONLY : RFLO_InitFlowSolver, RFLO_FlowSolverDummy, &
                            RFLO_FlowSolver, RFLO_UpdateInbuffGm, &
                            Fluid_finalize, Fluid_preHdfOutput, &
                            Fluid_postHdfOutput
  IMPLICIT NONE

  INCLUDE 'comf90.h'

  INTERFACE
    SUBROUTINE COM_set_pointer( attr,ptr,asso )
      USE ModRocstar, ONLY : t_globalGenx
      CHARACTER(*), INTENT(IN) :: attr
      TYPE(t_globalGenx), POINTER  :: ptr
      EXTERNAL asso
    END SUBROUTINE COM_set_pointer
  END INTERFACE

! ... parameters
  CHARACTER(*), INTENT(in) :: winName

! ... local variables
  LOGICAL :: isDummy

  INTEGER :: types(7)

  TYPE(t_globalGenx), POINTER  :: glb

!******************************************************************************

  ALLOCATE( glb)
  ALLOCATE( glb%global )

  isDummy            = TRIM(winName) == "RocfloDummy"
  glb%isDummy        = isDummy
  glb%global%winName = winName

! init fluids window

  CALL COM_new_window( winName )

! create an dataitem for global data

  CALL COM_new_dataitem( winName//'.global','w',COM_F90POINTER,1,'' )

  CALL COM_allocate_array( winName//'.global')

! solver initialization

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_MPI_COMM
  types(4) = COM_INTEGER
  types(5) = COM_STRING
  types(6) = COM_STRING
  types(7) = COM_INTEGER

  CALL COM_set_member_function( winName//'.initialize', &
          RFLO_InitFlowSolver,winName//'.global','biiiiii',types )

! solution update

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_DOUBLE_PRECISION
  types(4) = COM_INTEGER
  types(5) = COM_INTEGER

  IF (isDummy) THEN
    CALL COM_set_member_function( winName//'.update_solution', &
            RFLO_FlowSolverDummy,winName//'.global','biiii',types )
  ELSE
    CALL COM_set_member_function( winName//'.update_solution', &
            RFLO_FlowSolver,winName//'.global','biiii',types )
  ENDIF

! standalone mode

  CALL COM_set_member_function( winName//'.update_inbuff_gm_fluid', &
          RFLO_UpdateInbuffGm,winName//'.global','bi',types )

! shutdown

  CALL COM_set_member_function( winName//'.finalize', &
          Fluid_finalize,winName//'.global','b',types )

! IO

  CALL COM_set_member_function( winName//'.pre_hdf_output', &
          Fluid_preHdfOutput,winName//'.global','b',types )

  CALL COM_set_member_function( winName//'.post_hdf_output', &
          Fluid_postHdfOutput,winName//'.global','b',types )

! finish

  CALL COM_window_init_done( winName )

  CALL COM_set_pointer( WinName//'.global',glb,associate_pointer )

END SUBROUTINE Rocflo_load_module

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Rocflo_load_module.F90,v $
! Revision 1.6  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/12/08 19:55:46  wasistho
! added postHdfOutput
!
! Revision 1.3  2005/12/08 00:19:06  wasistho
! stored actual time averaged vars in hdf
!
! Revision 1.2  2005/02/03 23:04:38  jiao
! Changed the third argument to COM_MPI_COMM.
!
! Revision 1.1  2004/12/01 21:23:56  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/06/29 23:53:34  wasistho
! migrated to Roccom-3
!
! Revision 1.6  2003/12/07 04:06:26  jiao
! Changed the prototype of COM_set_pointer to be the same as COM_get_pointer,
! so that it takes ASSOCIATE_POINTER as the third argument. This change is
! needed to support PGI compilers.
!
! Revision 1.5  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/10/18 16:49:19  jblazek
! Changed parameter lists to some GenX routines.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







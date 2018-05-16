!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************

!!****
!!
!!  NAME
!!     implicit_finalize
!!
!!  FUNCTION
!!     Finalizes BlockSolve95 and deallocates implicit variables
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE implicit_finalize(global )

  USE implicit_global

  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: global
  INTEGER :: BS95debug

  IF ( global%debug_state ) THEN
     BS95debug = 1
  ELSE
     BS95debug = 0
  END IF

  CALL BS95FREE(BS95debug)
  CALL BS95FINALIZE(BS95debug)

  DEALLOCATE(cval_m)
  DEALLOCATE(aval_m)
  DEALLOCATE(rp_m)
  DEALLOCATE(cval_k)
  DEALLOCATE(aval_k)
  DEALLOCATE(rp_k)
  DEALLOCATE(cval_meff)
  DEALLOCATE(aval_meff)
  DEALLOCATE(rp_meff)
  DEALLOCATE(node_flag)
  DEALLOCATE(boundary_value)
  DEALLOCATE(Local2Global)
  DEALLOCATE(Global2Local)
  DEALLOCATE(NodeProc)
  DEALLOCATE(fext_imp)
  IF (nprocs > 1) THEN
     DEALLOCATE(CommProcs1)
     DEALLOCATE(NumCommNodes1)
     DEALLOCATE(CommNodes1)
     DEALLOCATE(CommProcsFrom1)
     DEALLOCATE(NumCommNodesFrom1)
     DEALLOCATE(CommProcs2)
     DEALLOCATE(NumCommNodes2)
     DEALLOCATE(CommNodes2)
     DEALLOCATE(CommProcsFrom2)
     DEALLOCATE(NumCommNodesFrom2)
  ENDIF


END SUBROUTINE implicit_finalize


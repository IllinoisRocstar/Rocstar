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
SUBROUTINE Rocfrac_load_module( module_name)
  USE RocFracMain
  USE RocFracSubInterface
  USE IntegralCheck
  IMPLICIT NONE
  INCLUDE 'comf90.h'
  INCLUDE 'rocmanf90.h'
    
  CHARACTER(*), INTENT(IN) :: module_name

  LOGICAL   :: isDummy
  INTEGER :: types(7)
  TYPE(ROCFRAC_GLOBAL), POINTER :: glb
    
  INTERFACE 
     SUBROUTINE COM_set_pointer( attr, ptr, asso)
       USE      ROCSTAR_RocFrac
       CHARACTER(*), INTENT(IN) :: attr
       TYPE(ROCFRAC_GLOBAL), POINTER :: ptr
       EXTERNAL asso
     END SUBROUTINE COM_set_pointer
  END INTERFACE
    
  isDummy = TRIM(module_name) == "RocfracDummy"
  
  ALLOCATE( glb)
  ! dummy rocfrac flag
  glb%iDummyRocfrac = isDummy
    
  CALL COM_new_window( module_name)
  
  CALL COM_new_dataitem( module_name//".global", 'w', COM_F90POINTER, 1, '')
  CALL COM_resize_array( module_name//".global")
  
  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_MPI_COMM
  types(4) = COM_INTEGER 
  types(5) = COM_STRING
  types(6) = COM_STRING
  types(7) = COM_INTEGER
 
  CALL COM_set_member_function( module_name//".initialize", &
       RocFracInitialize, module_name//".global", "biiiiii", types)

  CALL COM_set_member_function( module_name//".finalize", &
       RocFracFinalize, module_name//".global", "b", types)

  types(3) = COM_DOUBLE_PRECISION
  types(4) = COM_INTEGER
  IF ( .NOT. isDummy) THEN
     CALL COM_set_member_function( module_name//".update_solution", &
          RocFracSoln, module_name//".global", "biii", types)
  ELSE
     CALL COM_set_member_function( module_name//".update_solution", &
          RocFracInterfaceUpdate, module_name//".global", "biii", types)
  END IF

!!! This is for running Rocfrac in standalone mode
  CALL COM_set_member_function( module_name//".update_inbuff_bc_solid", &
       RocFracUpdateInbuff, module_name//".global", "bi", types)

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION

  CALL COM_set_member_function( module_name//'.'//MAN_COMP_INTEG_NAME, &
       CheckIntegral, module_name//".global", "bo", types)
  CALL COM_window_init_done( module_name)

  CALL COM_set_pointer( module_name//".global", glb, associate_pointer)

END SUBROUTINE Rocfrac_load_module

SUBROUTINE Rocfrac_unload_module( module_name) 
  USE ROCSTAR_RocFrac
  IMPLICIT NONE
  INCLUDE 'comf90.h'

  CHARACTER(*), INTENT(IN) :: module_name

  TYPE(ROCFRAC_GLOBAL), POINTER :: glb

  INTERFACE 
     SUBROUTINE COM_get_pointer( attr, ptr, asso)
       USE RocFracMain
       CHARACTER(*), INTENT(IN) :: attr
       TYPE(ROCFRAC_GLOBAL), POINTER :: ptr
       EXTERNAL asso
     END SUBROUTINE COM_get_pointer
  END INTERFACE

  CALL COM_get_pointer( module_name//".global", glb, associate_pointer)
  DEALLOCATE( glb)

  CALL COM_delete_window( module_name)

END SUBROUTINE Rocfrac_unload_module


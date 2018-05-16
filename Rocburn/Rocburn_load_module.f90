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
!
! ---------------------------------------------------------------------------
!
!  Purpose: Registering subroutines of ROCBURN with Roccom. This (hopefully)
!           is the only file need to be changed to add another 1D model.
!
!  Author:            X. Jiao
!
!  Creation Date:     Aug. 30, 2002
!
!  Modifications:
!
!    No.     Date         Programmer    Description
!
! ---------------------------------------------------------------------------
!

MODULE M_INIT_1DMODULES
  USE M_ROCBURN_2D
  IMPLICIT NONE

  INTERFACE 
     SUBROUTINE COM_set_data( attr, ptr)  ! Must be explicitly declared here
       USE M_ROCBURN_INTERFACE_DATA
       CHARACTER(*), INTENT(IN)   :: attr
       EXTERNAL  ptr
     END SUBROUTINE COM_set_data
  END INTERFACE

CONTAINS  

!
! The following subroutines registers subroutines of indiviudal 1-D burning 
!   rate model to Roccom. Rocman calls corresponding subroutine
!   related to the burning rate model used.
!

  SUBROUTINE ROCBURN_INIT_FUNCS_APN( mname) 
    USE M_ROCBURN_1D_APN
    CHARACTER(*), INTENT(IN) :: mname
    
    CALL COM_set_external( mname//".init_0d", 0, INITIALIZE_0D)
    CALL COM_set_external( mname//".init_1d", 0, INITIALIZE_1D)
    CALL COM_set_external( mname//".finalize_0d", 0, FINALIZE_0D)
    CALL COM_set_external( mname//".get_burn_rate", 0, GET_BURNING_RATE_1D)

  END SUBROUTINE ROCBURN_INIT_FUNCS_APN

  SUBROUTINE ROCBURN_INIT_FUNCS_PY( mname) 
    USE M_ROCBURN_1D_PY
    CHARACTER(*), INTENT(IN) :: mname

    CALL COM_set_external( mname//".init_0d", 0, INITIALIZE_0D)
    CALL COM_set_external( mname//".init_1d", 0, INITIALIZE_1D)
    CALL COM_set_external( mname//".finalize_0d", 0, FINALIZE_0D)
    CALL COM_set_external( mname//".get_film_coeff", 0, GET_FILM_COEFF_1D)
    CALL COM_set_external( mname//".get_time_step", 0, GET_TIME_STEP_1D)
    CALL COM_set_external( mname//".get_burn_rate", 0, GET_BURNING_RATE_1D)

  END SUBROUTINE ROCBURN_INIT_FUNCS_PY

  SUBROUTINE ROCBURN_INIT_FUNCS_ZN( mname) 
    USE M_ROCBURN_1D_ZN
    CHARACTER(*), INTENT(IN) :: mname

    CALL COM_set_external( mname//".init_0d", 0, INITIALIZE_0D)
    CALL COM_set_external( mname//".init_1d", 0, INITIALIZE_1D)
    CALL COM_set_external( mname//".finalize_0d", 0, FINALIZE_0D)
    CALL COM_set_external( mname//".get_film_coeff", 0, GET_FILM_COEFF_1D)
    CALL COM_set_external( mname//".get_time_step", 0, GET_TIME_STEP_1D)
    CALL COM_set_external( mname//".get_burn_rate", 0, GET_BURNING_RATE_1D)

  END SUBROUTINE ROCBURN_INIT_FUNCS_ZN
END MODULE M_INIT_1DMODULES

SUBROUTINE ROCBURN_LOAD_MODULE( module_name)

  USE M_INIT_1DMODULES
  USE M_ROCBURN_INTERFACE_DATA, ONLY : associate_pointer
  IMPLICIT NONE
 
  CHARACTER(*), INTENT(IN)  :: module_name
  TYPE(list_block), POINTER :: G_b
  INTEGER                   :: types(7)

  INTERFACE 
     SUBROUTINE COM_set_pointer( attr, ptr, asso)
       USE       M_ROCBURN_INTERFACE_DATA, ONLY : list_block
       CHARACTER(*), INTENT(IN) :: attr
       TYPE(list_block), POINTER :: ptr
       EXTERNAL asso
     END SUBROUTINE COM_set_pointer
  END INTERFACE

  ALLOCATE( G_b)
  G_b%mname = module_name
  G_b%TBL_flag = NO_TBL
  IF ( module_name == "RocburnAPN") THEN
     G_b%burn_model = MODEL_APN
  ELSE IF ( module_name == "RocburnPY") THEN
     G_b%burn_model = MODEL_PY
  ELSE IF ( module_name == "RocburnZN") THEN
     G_b%burn_model = MODEL_ZN
  ELSE
     PRINT *, "Rocburn-2D: Unknown module name", module_name
     PRINT *, "Rocburn-2D: Use APN instead", module_name
     G_b%burn_model = MODEL_APN
  END IF

  CALL COM_new_window( module_name)
!!! Create an dataitem for global data
  CALL COM_new_dataitem( module_name//".global", 'w', COM_F90POINTER, 1, '')
  CALL COM_resize_array( module_name//".global")

  CALL COM_new_dataitem( module_name//".init_0d", 'w', COM_VOID, 1, '')
  CALL COM_new_dataitem( module_name//".init_1d", 'w', COM_VOID, 1, '')
  CALL COM_new_dataitem( module_name//".finalize_0d", 'w', COM_VOID, 1, '')
  CALL COM_new_dataitem( module_name//".get_film_coeff", 'w', COM_VOID, 1, '')
  CALL COM_new_dataitem( module_name//".get_time_step", 'w', COM_VOID, 1, '')
  CALL COM_new_dataitem( module_name//".get_burn_rate", 'w', COM_VOID, 1, '')

  CALL COM_resize_array( module_name//".init_0d")
  CALL COM_resize_array( module_name//".init_1d")
  CALL COM_resize_array( module_name//".finalize_0d")
  CALL COM_resize_array( module_name//".get_film_coeff")
  CALL COM_resize_array( module_name//".get_time_step")
  CALL COM_resize_array( module_name//".get_burn_rate")


!!! Now initialize the 1D module
  IF ( G_b%burn_model == MODEL_PY) THEN
     CALL ROCBURN_INIT_FUNCS_PY( module_name)
  ELSE IF ( G_b%burn_model == MODEL_ZN) THEN
     CALL ROCBURN_INIT_FUNCS_ZN( module_name)
  ELSE 
     CALL ROCBURN_INIT_FUNCS_APN( module_name)
  END IF
  
  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_MPI_COMM
  types(4) = COM_INTEGER
  types(5) = COM_STRING
  types(6) = COM_STRING
  types(7) = COM_INTEGER

  CALL COM_set_member_function( module_name//".initialize", &
       INIT_WRAPPER, module_name//".global", "biiiiii", types)
  
  types(1) = COM_F90POINTER
  types(2) = COM_INTEGER
  types(3) = COM_STRING
  types(4) = COM_STRING
  types(5) = COM_RAWDATA
  types(6) = COM_RAWDATA
  types(7) = COM_INTEGER

  CALL COM_set_member_function( module_name//".init_internal", &
       INITIALIZE, module_name//".global", "biiiiii", types)
  
  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_DOUBLE_PRECISION
  types(4) = COM_INTEGER

  CALL COM_set_member_function( module_name//".update_solution",  &
       UPDATE_WRAPPER, module_name//".global", "biii", types)

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_DOUBLE_PRECISION
  types(4) = COM_INTEGER
  types(5) = COM_RAWDATA
  types(6) = COM_RAWDATA
  types(7) = COM_RAWDATA

  CALL COM_set_member_function( module_name//".update_internal",  &
       UPDATE, module_name//".global", "biiiiii", types)


  types(1) = COM_F90POINTER
  types(2) = COM_RAWDATA
  CALL COM_set_member_function( module_name//".finalize", &
       FINALIZE, module_name//".global", "b", types)

  CALL COM_window_init_done( module_name)

  G_b%INIT = COM_get_function_handle(module_name//".init_internal")
  G_b%UPDATE = COM_get_function_handle(module_name//".update_internal")
  G_b%INIT_1D = COM_get_dataitem_handle(module_name//".init_1d")
  G_b%INIT_0D = COM_get_dataitem_handle(module_name//".init_0d")
  G_b%GET_FILM_COEFF = COM_get_dataitem_handle(module_name//".get_film_coeff")
  G_b%GET_BURN_RATE  = COM_get_dataitem_handle(module_name//".get_burn_rate")
  G_b%GET_TIME_STEP = COM_get_dataitem_handle(module_name//".get_time_step")

  CALL COM_set_pointer( module_name//".global", G_b, associate_pointer)
END SUBROUTINE ROCBURN_LOAD_MODULE

SUBROUTINE ROCBURN_UNLOAD_MODULE( module_name)
  USE M_ROCBURN_INTERFACE_DATA, ONLY : list_block, associate_pointer
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: module_name

  TYPE( list_block), POINTER  :: glb
  
  INTERFACE 
     SUBROUTINE COM_get_pointer( attr, ptr, asso)
       USE M_ROCBURN_INTERFACE_DATA
       CHARACTER(*), INTENT(IN)  :: attr
       TYPE(list_block), POINTER :: ptr
       EXTERNAL asso
     END SUBROUTINE COM_get_pointer
  END INTERFACE

  CALL COM_get_pointer( module_name//".global", glb, associate_pointer)
  DEALLOCATE( glb)

  CALL COM_delete_window( module_name)

END SUBROUTINE ROCBURN_UNLOAD_MODULE







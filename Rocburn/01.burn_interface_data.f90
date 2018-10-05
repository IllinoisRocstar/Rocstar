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
MODULE M_ROCBURN_INTERFACE_DATA
!
!  This module contains the interface variables needed for all 
!  burning rate models.  
!
!  The burning rate models available are:
!
!    Rocburn_APN  : Quasi-steady model using rb=a*P^n
!    Rcoburn_PY   : Model implemented by  Luca Massa
!                   uses a pyrolysis law and large activation energy (Buckmaster)
!    Rocburn_ZN   : Model implemented by  Tang and Brewster 
!                         using Zeldovich-Novozhilov approach
!
!  Author:          K-C Tang, L. Massa, X. Jiao
!

  IMPLICIT NONE


!  ===================
!  D A T A   T Y P E S
!  ===================
  INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(P=14,R=30)

!
!  ---------------
!  Block data type
!  ---------------

  TYPE, PUBLIC :: block

!
!           data for initialization
!
     
     INTEGER                    :: iblock           ! block id
     INTEGER                    :: nfaces           ! number of faces 

!
!  ---------------------------------------------
!       incoming data from ROCMAN
!
     REAL(DBL), POINTER         :: coor(:,:)
     REAL(DBL), POINTER         :: pres(:)
     REAL(DBL), POINTER         :: qr(:), qc(:)     ! not used for APN
     REAL(DBL), POINTER         :: rhoc(:)          ! not used for APN   
!    REAL(DBL), POINTER         :: To(:)            ! from solid if available

!       incoming data for specific burning rate model

     REAL(DBL), POINTER         :: Tg(:)            ! Rocburn_PY
     INTEGER, POINTER           :: burn_flag(:)     ! Rocburn_PY

!
!
!  ---------------------------------------------
!       outgoing data to ROCMAN
!
     REAL(DBL), POINTER         :: rb(:)
     REAL(DBL), POINTER         :: Tf(:)

!  ---------------------------------------------
!           
!       internal data storage for burning rate models 
!         
!
!          independent variables; old state variables
!
     REAL(DBL), POINTER         :: qc_old(:), qr_old(:) 
     REAL(DBL), POINTER         :: pres_old(:), Tg_old(:)
     REAL(DBL), POINTER         :: rhoc_old(:)
!    REAL(DBL), POINTER         :: To_old(:)  

!
!          dependent variables
!

     REAL(DBL), POINTER         :: temp(:,:)
!
!          data storage for specific burning rate model --- dependent varialbes
!


     REAL(DBL), POINTER         :: Toa(:)           ! Rocburn_ZN
     REAL(DBL), POINTER         :: fr(:)            ! Rocburn_ZN
     REAL(DBL), POINTER         :: dist(:)          ! Rocburn_PY

  END TYPE block

! This is a placeholder for the type defined in 1D modules
  TYPE, PUBLIC :: G_BURN_1D
     INTEGER :: buf(4096)
  END TYPE G_BURN_1D

!
!  -----------------------
!  L I N K E D   L I S T S
!  -----------------------

  TYPE, PUBLIC :: list_block
     TYPE(block), POINTER      :: blocks(:)

     INTEGER                   :: MPI_COMM_ROCBURN, rank
     INTEGER                   :: burn_model, TBL_flag, burn_iter, burn_cell, total_cell

!           function handles introduced for Roccom 3
     INTEGER                    :: INIT, UPDATE
     INTEGER                    :: INIT_0D, INIT_1D
     INTEGER                    :: GET_TIME_STEP, GET_BURN_RATE, GET_FILM_COEFF

     REAL(DBL)                 :: To_read, pseudo_time
     REAL(DBL), POINTER        :: Tn(:)            ! Buffer for 1D rocburn

     CHARACTER(LEN=80)         :: mname

     TYPE(G_BURN_1D), POINTER  :: g_1d
  END TYPE list_block

CONTAINS

  SUBROUTINE associate_pointer( attr, ptr)
    TYPE(list_block), POINTER :: attr, ptr
    ptr => attr
  END SUBROUTINE ASSOCIATE_POINTER

END MODULE M_ROCBURN_INTERFACE_DATA







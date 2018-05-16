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
! Purpose: Compute information for proper use of DCUHRE
!
! Description: None.
!
! Input: 
!   NDIM        Number of dimensions
!   NF          Number of functions
!   KEY         Type of integration rule
!   MAXCLS      Maximum number of calls allowed (NOTE also output argument)
!
! Output:
!   MAXCLS      Maximum number of calls allowed (NOTE also input argument)
!   NW          Size of work array
!
! Notes:
!
!******************************************************************************
!
! $Id: RFLU_ComputeDCUHREInfo.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

  SUBROUTINE RFLU_ComputeDCUHREInfo(global,NDIM,NF,KEY,MAXCLS,NW)

    USE ModGlobal, ONLY: t_global
    USE ModError

    IMPLICIT NONE

! - parameters      

    INTEGER, INTENT(IN) :: KEY,NDIM,NF
    INTEGER, INTENT(OUT) :: NW
    INTEGER, INTENT(INOUT) :: MAXCLS    
    TYPE(t_global), POINTER :: global

! - locals

    INTEGER :: MAXSUB,NUM

!******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeDCUHREInfo',&
  'RFLU_ComputeDCUHREInfo.F90')

! - Compute information ------------------------------------------------------

    SELECT CASE (KEY) ! Select according to integration rule
      CASE (0)  
        NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) &
            + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
      CASE (4)
        NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! KEY

    MAXCLS = MAX(MAXCLS,4*NUM) ! documentation recommends MAXCLS >= 3*NUM      

    MAXSUB = (MAXCLS - NUM)/(2*NUM) + 1
    NW     = MAXSUB*(2*NDIM+2*NF+2) + 17*NF + 1

    CALL DeregisterFunction(global)

!******************************************************************************

  END SUBROUTINE RFLU_ComputeDCUHREInfo
  
!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeDCUHREInfo.F90,v $
! Revision 1.4  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2002/09/09 14:15:01  haselbac
! global now under regions
!
! Revision 1.1  2002/07/25 14:34:59  haselbac
! Initial revision
!
!******************************************************************************
  








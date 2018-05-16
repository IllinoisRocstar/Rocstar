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
! Purpose: register portable random number variables with GenX.
!
! Description: none.
!
! Input: regions = patches and region (volume) variables
!
! Output: to Roccom via randInitGenxInterface.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RandInitGenxInterface.F90,v 1.3 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE randInitGenxInterface( regions, wins, winv )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
#endif
  USE ModError
  USE ModParameters
  USE ModRandom
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  CHARACTER(CHRLEN) :: wins, winv
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER :: pid 
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'randInitGenxInterface',&
  'RandInitRocstarInterface.F90' )

! input data (currently none) -------------------------------------------------

! output surface data (currently none) ----------------------------------------
 
! output restart data ---------------------------------------------------------

  CALL COM_new_dataitem( TRIM(winv)//'.rand', 'p',         &
                          COM_INTEGER,1, '')
  CALL COM_set_size( TRIM(winv)//'.rand', 0, RAND_TOTAL_SIZE)

! store pointers to variables, loop over all regions --------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE          ) THEN    ! on my processor

! --- volume data -------------------------------------------------------------

      pid  = iReg*REGOFF
      CALL COM_set_array( TRIM(winv)//'.rand',pid,regions(iReg)%randData%mt(0))

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE randInitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RandInitGenxInterface.F90,v $
! Revision 1.3  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:47  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/06/30 04:07:48  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.4  2004/06/29 23:53:25  wasistho
! migrated to Roccom-3
!
! Revision 1.3  2004/03/01 23:47:35  jiao
! Changed the F90 implementation for COM_init_dataitem and COM_init_mesh
! not to require registered scalars to be F90 pointers.
!
! Revision 1.2  2004/02/20 00:49:04  jiao
! Changed to use COM_INIT_ATTRIBUTE_SPECIAL in order to compile with older IBM compilers.
!
! Revision 1.1  2003/11/21 22:20:04  fnajjar
! Initial import
!
!******************************************************************************








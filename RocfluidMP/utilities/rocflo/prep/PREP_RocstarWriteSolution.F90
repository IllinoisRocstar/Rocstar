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
! Purpose: write flow solution to file for one region.
!
! Description: solution is written in user selected Genx format
!
! Input: gridLevel  = initial grid level
!        region     = dimensions and cons. variables of current region
!        iReg       = region number
!        wins, winv = surface and volume window names
!
! Output: to file through Rocout.
!
! Notes: solution is stored only for the current grid level; it is also
!        stored for all dummy cells. All regions are written into one file.
!        There is no transfer of data from other processors.
!
!******************************************************************************
!
! $Id: PREP_GenxWriteSolution.F90,v 1.4 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GenxWriteSolution( gridLevel,iReg,region,wins,winv )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
#ifdef GENX
  USE PREP_ModInterfaces, ONLY : GenxWriteRocinout
#endif
  USE ModParameters
  IMPLICIT NONE
  INCLUDE "comf90.h"

! ... parameters
  INTEGER :: gridLevel, iReg
  TYPE(t_region)    :: region
  CHARACTER(CHRLEN) :: wins, winv
  CHARACTER(CHRLEN) :: fname

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  TYPE(t_global), POINTER :: global
  INTEGER :: write_attr, set_option, surf_all, vol_all, errFlg

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'GenxWriteSolution',&
  'PREP_RocstarWriteSolution.F90' )

! obtain function handle ------------------------------------------------------

  write_attr = COM_get_function_handle( 'OUT.write_dataitem')
  set_option = COM_get_function_handle( 'OUT.set_option')

  IF ( iReg==1) THEN
     CALL GenxWriteRocinout( global )
  ENDIF
  CALL COM_call_function( set_option, 2, 'mode', 'w')

! do not append process rank --------------------------------------------------

  CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write volume window ---------------------------------------------------------

  vol_all = Com_get_dataitem_handle( TRIM(winv)//'.all')

  WRITE(fname,'(A,I5.5)') '../Rocin/fluid_',iReg 
  CALL COM_call_function( write_attr, 4, TRIM(fname), vol_all, &
       "fluid", "00.000000")

! write surface window --------------------------------------------------------

  surf_all = Com_get_dataitem_handle( TRIM(wins)//'.all')

  WRITE(fname,'(A,I5.5)') '../Rocin/ifluid_',iReg
  CALL COM_call_function( write_attr, 4, TRIM(fname), surf_all, &
       "ifluid", "00.000000")

! delete volume and surface windows -------------------------------------------

  CALL COM_delete_window( TRIM(winv))
  CALL COM_delete_window( TRIM(wins))

! deallocate arrays -----------------------------------------------------------

  DEALLOCATE( region%levels(gridLevel)%grid%xyz,stat=errFlg )
  global%error = errFlg
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE( region%levels(gridLevel)%mixt%cv,stat=errFlg )
  global%error = errFlg
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DO iPatch=1,region%nPatches
    DEALLOCATE( region%levels(gridLevel)%patches(iPatch)%surfCoord,stat=errFlg )
    global%error = errFlg
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

    DEALLOCATE( region%levels(gridLevel)%patches(iPatch)%bcFlag,stat=errFlg )
    global%error = errFlg
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GenxWriteSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_GenxWriteSolution.F90,v $
! Revision 1.4  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:28:46  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.5  2004/10/12 04:32:52  wasistho
! split to one block per HDF file
!
! Revision 1.4  2004/07/27 03:33:23  wasistho
! restructured directories
!
! Revision 1.3  2004/07/23 23:24:11  wasistho
! moved creation of Rocin under Rocflo
!
! Revision 1.2  2004/07/23 04:31:03  wasistho
! Genx: readin from Rocin, standalone: read .inp file i.o. command line input
!
! Revision 1.1  2004/06/30 00:06:05  wasistho
! initial import for GEN3
!
!
!******************************************************************************








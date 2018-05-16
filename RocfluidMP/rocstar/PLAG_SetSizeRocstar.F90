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
! Purpose: Set size of arrays in Genx.
!
! Description: none.
!
! Input: regions    = data of all regions,
!        wins, winp = GenX window registrations.
!
! Output: to Roccom.
!
! Notes: Surface registration for Tiles works only for External coupled bc.
!        Need to activate for both Internal and External bc.
!
!******************************************************************************
!
! $Id: PLAG_SetSizeGenx.F90,v 1.5 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_SetSizeGenx( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE ModPartLag, ONLY    : t_plag
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters

  TYPE(t_region) :: region

! ... loop variables

! ... local variables

  INTEGER :: iLev, pid
  TYPE(t_global),     POINTER :: pGlobal
  TYPE(t_plag)  ,     POINTER :: pPlag

!******************************************************************************

  pGlobal => region%global

  CALL RegisterFunction( pGlobal,'PLAG_SetSizeGenx',&
  'PLAG_SetSizeRocstar.F90' )

!------------------------------------------------------------------------------
! Loop over all regions 
!------------------------------------------------------------------------------

  iLev   = region%currLevel
  pPlag => region%levels(iLev)%plag
  pid = region%iRegionGlobal*REGOFF+1

!------------------------------------------------------------------------------
! COM_set_size procedure must be called in RFLO_sendBoundaryValues
! if nPcls has changed.
!------------------------------------------------------------------------------
      
  CALL COM_set_size( TRIM(pGlobal%winp)//'.nc',pid,pPlag%nPcls)

!------------------------------------------------------------------------------
! finalize 
!------------------------------------------------------------------------------

  CALL DeregisterFunction( pGlobal )

END SUBROUTINE PLAG_SetSizeGenx

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_SetSizeGenx.F90,v $
! Revision 1.5  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/11/16 08:51:13  jiao
! Fixed pane ID.
!
! Revision 1.2  2004/07/02 22:47:43  jiao
! Added definitions of iLev and pid.
!
! Revision 1.1  2004/07/02 22:08:04  fnajjar
! Initial import for Roccom3
!
!******************************************************************************








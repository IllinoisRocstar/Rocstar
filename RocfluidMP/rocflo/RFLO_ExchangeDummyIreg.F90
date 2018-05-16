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
! Purpose: exchange values between interior and dummy cells of the
!          corresponding patches of the adjacent regions (both regions
!          are on the same processor).
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: regions%levels%mixt%cv = flow variables in dummy cells.
!
! Notes: intended for non-conforming region interface with completely
!        unrelated grids.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeDummyIreg.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeDummyIreg( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : 
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables


! ... local variables

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ExchangeDummyIreg',&
  'RFLO_ExchangeDummyIreg.F90' )



  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ExchangeDummyIreg

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeDummyIreg.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
!******************************************************************************








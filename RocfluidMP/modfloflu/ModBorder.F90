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
! ******************************************************************************
!
! Purpose: Definition of derived data type for borders.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModBorder.F90,v 1.9 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModBorder

  USE ModDataTypes
  
  IMPLICIT NONE

! ******************************************************************************
! Type definition
! ******************************************************************************

  TYPE t_border_data
    INTEGER :: sendRequest,sendRequestCount,sendRequestInt,tag,tagCount,tagInt
    INTEGER, DIMENSION(:,:), POINTER :: recvBuffInt,sendBuffInt
    REAL(RFREAL), DIMENSION(:,:), POINTER :: recvBuff,sendBuff
  END TYPE t_border_data

  TYPE t_border
    INTEGER :: iBorder,iProc,iRegionGlobal,iRegionLocal
    INTEGER :: nCellsRecv,nCellsSend,nVertRecv,nVertSend,nVertShared
#ifdef PLAG
    INTEGER :: nPclsRecv,nPclsSend,nPclsSendMax
#endif
    INTEGER, DIMENSION(:), POINTER :: icgRecv,icgSend,ivgRecv,ivgSend,ivgShared
#ifdef PLAG
    INTEGER, DIMENSION(:,:), POINTER :: iPclSend
#endif
    TYPE(t_border_data) :: mixt,spec,plag
  END TYPE t_border

END MODULE ModBorder

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModBorder.F90,v $
! Revision 1.9  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2005/12/14 21:19:16  fnajjar
! Added nPclsSendMax for dynamic iPclSend
!
! Revision 1.6  2005/12/13 23:06:16  fnajjar
! Added tags and sendRequests pertinent for PLAG
!
! Revision 1.5  2005/05/18 22:05:26  fnajjar
! Added integer vars, fixed bug in declaration
!
! Revision 1.4  2005/04/30 13:47:51  haselbac
! Added arrays for parallelization of particle module
!
! Revision 1.3  2005/04/15 15:06:26  haselbac
! Added data arrays and variables
!
! Revision 1.2  2005/01/14 21:14:20  haselbac
! Added iProc
!
! Revision 1.1  2004/12/04 03:43:41  haselbac
! Initial revision
!
! ******************************************************************************







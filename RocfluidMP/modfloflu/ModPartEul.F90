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
! Purpose: define data types related to Eulerian particles.
!
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModPartEul.F90,v 1.12 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModPartEul

  USE ModDataTypes
  USE ModMaterials, ONLY : t_material
  IMPLICIT NONE

! input -----------------------------------------------------------------------

  TYPE t_peul_input
    INTEGER :: nPtypes, readStatus
    REAL(RFREAL) :: smoocf
    LOGICAL :: constInit ! if allowed to initialize to constant for t > 0
    TYPE(t_peul_ptype), POINTER :: ptypes(:)
  END TYPE t_peul_input

! particle type ---------------------------------------------------------------

  TYPE t_peul_ptype
! input variables
    TYPE(t_material), POINTER :: material
    REAL(RFREAL) :: diam,puff,initc
    REAL(RFREAL) :: Sc,vis2,vis4
    INTEGER      :: negReport,clipModel,methodV
! derived variables
    REAL(RFREAL) :: denseff,voleff,volmat,tauVcoef,maxConc
  END TYPE t_peul_ptype

! data ------------------------------------------------------------------------

  TYPE t_peul
    INTEGER :: nCv1, nDv1, nTv1 ! number of variables per particle type
    INTEGER :: nCv,  nDv,  nTv  ! total number of variables
    INTEGER :: nRequests

    INTEGER,      POINTER :: requests(:)
    REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:), dv(:,:), tv(:,:)
    REAL(RFREAL), POINTER :: rhs(:,:), rhsSum(:,:), diss(:,:), fterm(:,:)

#ifdef RFLO
    REAL(RFREAL), POINTER :: srad(:,:), epsIrs(:,:)
#endif
  END TYPE t_peul

! communication data structure -------------------------------------------------

  TYPE t_buffer_peul
    INTEGER :: iRequest, nRecvBuff, nSendBuff
    REAL(RFREAL), POINTER, DIMENSION(:) :: recvBuff, sendBuff
  END TYPE t_buffer_peul

END MODULE ModPartEul

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModPartEul.F90,v $
! Revision 1.12  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2004/05/03 15:09:41  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.9  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.8  2004/03/02 21:44:51  jferry
! Added clipping options
!
! Revision 1.7  2003/04/14 16:33:52  jferry
! added option to initialize to constant for t > 0
!
! Revision 1.6  2003/04/09 14:08:13  fnajjar
! Added datastructure for MPI communication
!
! Revision 1.5  2003/04/07 18:25:52  jferry
! added initial concentration for smoke to data structure
!
! Revision 1.4  2003/03/11 16:01:23  jferry
! Created data type for material properties
!
! Revision 1.3  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************







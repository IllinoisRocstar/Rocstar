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
! Purpose: encapsulate MPI header file, define master process.
!
! Description: none
!
! Notes: used in order to avoid the INCLUDE statement in the source code
!        (not standard in Fortran 90)
!
! ******************************************************************************
!
! $Id: ModMPI.F90,v 1.16 2009/09/07 13:10:54 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModMPI

! Trying to fix ill-conceived tag management.  Note that the PLAG_TAG_SHIFT 
! exceeds the maximum value tags are allowed to take in the MPI standard.  Tags
! only need to be unique among pairs of processors.  All RocfluidMP comms are
! pair-wise. -mtc
!
! 
!  INTEGER, PARAMETER :: MASTERPROC   = 0, &   ! master process
!                        MPI_PATCHOFF = 100    ! offset for patch numbers (tag)
!
!  INTEGER, PARAMETER :: BASE_TAG_SHIFT = 10000, & ! MPI tag shifts 
!                        PEUL_TAG_SHIFT = 10001, & ! all must be positive 
!                        TURB_TAG_SHIFT = 20001, & ! shift = n*base+1, n=1,2,..
!                        RADI_TAG_SHIFT = 30001, & ! 
!                        PLAG_TAG_SHIFT = 40001    ! PLAG always last
!  
! 
  INTEGER, PARAMETER :: MASTERPROC   = 0, &    ! master process
                        MPI_PATCHOFF = 100       ! offset for patch numbers (tag)

  INTEGER, PARAMETER :: BASE_TAG_SHIFT = 512, & ! MPI tag shifts 
                        PEUL_TAG_SHIFT = 513, & ! all must be positive 
                        TURB_TAG_SHIFT = 2049, & ! shift = n*base+1, n=1,2,..
                        RADI_TAG_SHIFT = 4097, & ! 
                        PLAG_TAG_SHIFT = 8193    ! PLAG always last

#ifdef RFLO
#ifdef MPI
  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: MPI_RFREAL = MPI_DOUBLE_PRECISION  ! goes with RFREAL
#endif
#endif
#ifdef RFLU
  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: MPI_RFREAL = MPI_DOUBLE_PRECISION  ! goes with RFREAL
#endif

END MODULE ModMPI

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModMPI.F90,v $
! Revision 1.16  2009/09/07 13:10:54  mtcampbe
! Increased tag shift for patches, should be OK now that values will be
! wrapped around above 32768.
!
! Revision 1.15  2009/08/12 04:15:57  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.14  2009/04/10 14:21:10  mtcampbe
! Updated with working values for tag/shifts.
!
! Revision 1.13  2009/04/07 15:11:04  mtcampbe
! Trying to fix ill-conceived tag management strategy
!
! Revision 1.12  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2005/04/15 15:34:09  haselbac
! Bug fix
!
! Revision 1.9  2005/04/15 15:30:28  haselbac
! Changed so that mpif.h always included for RFLU
!
! Revision 1.8  2004/10/19 19:28:58  haselbac
! Cosmetics only
!
! Revision 1.7  2004/09/27 22:47:00  wasistho
! added RADI_TAG_SHIFT
!
! Revision 1.6  2004/03/06 02:33:04  wasistho
! moved mpi tag shifts from ModParameters to ModMPI
!
! Revision 1.5  2002/03/21 18:07:15  jblazek
! Added check of MPI_PATCHOFF (for tags).
!
! Revision 1.4  2002/03/18 23:07:19  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/01/31 20:23:59  jblazek
! Added treatment of edge & corner cells.
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************







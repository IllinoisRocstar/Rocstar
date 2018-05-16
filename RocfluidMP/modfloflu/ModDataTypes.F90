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
! Purpose: Set the precision of reals and the length of character variables.
!
! Description: None.
!
! Notes: 
!   1. When changing RFREAL do not forget to adjust MPI_RFREAL in ModMPI.F90.
!
! ******************************************************************************
!
! $Id: ModDataTypes.F90,v 1.5 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModDataTypes
  INTEGER, PARAMETER :: SPREAL = SELECTED_REAL_KIND( 6, 37), & 
                        RFREAL = SELECTED_REAL_KIND(15,307), & 
                        CHRLEN = 80
END MODULE ModDataTypes

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModDataTypes.F90,v $
! Revision 1.5  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/10/05 20:05:56  haselbac
! Added SPREAL for ENSIGHT filter
!
! Revision 1.2  2002/03/18 23:07:19  jblazek
! Finished multiblock and MPI.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************







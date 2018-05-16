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
! Purpose: map regions to processors; assign local region numbers
!          on each processor.
!
! Description: none.
!
! Input: global = no. of processors, mapping type.
!
! Output: regions = processor ID, local region number.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_MapRegionsProcessors.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_MapRegionsProcessors( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModParameters
  USE ModError
  USE ModMPI
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iproc, iReg

! ... local variables
  INTEGER :: iRegBeg, iRegEnd, local, iprat, ipmod

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MapRegionsProcessors',&
  'RFLO_MapRegionsProcessors.F90' )

! "manual" region to processor mapping

  IF (global%nRegionsProc < 0) CALL ErrorStop( global,ERR_NO_PROCMAP,__LINE__ )

! automatic region to processor mapping

  IF (global%nRegionsProc == 0) THEN

    iprat = MAX(1,global%nRegions/global%nProcAlloc)
    ipmod = MOD(global%nRegions,global%nProcAlloc)

    iRegBeg = 1
    iRegEnd = iprat
    IF (ipmod > 0) THEN
      iRegEnd = iRegEnd + 1
      ipmod   = ipmod - 1
    ENDIF
    DO iproc=0,global%nProcAlloc-1
      local = 0
      DO iReg=iRegBeg,iRegEnd      ! regions on iproc
        local = local + 1
        IF (global%nProcAlloc>1 .AND. local>=MPI_PATCHOFF) &
          CALL ErrorStop( global,ERR_PATCH_OFFSET,__LINE__ )
        IF (global%nProcAlloc > 1) THEN
          regions(iReg)%procid = iproc
        ELSE
          regions(iReg)%procid = MASTERPROC
        ENDIF
        regions(iReg)%localNumber = local
      ENDDO
      iRegBeg = iRegEnd + 1
      iRegEnd = iRegBeg + iprat - 1
      IF (ipmod > 0) THEN
        iRegEnd = iRegEnd + 1
        ipmod   = ipmod - 1
      ENDIF
    ENDDO

  ENDIF

! specified # of regions per processor

  IF (global%nRegionsProc > 0) THEN

    IF (global%nProcAlloc*global%nRegionsProc /= global%nRegions) &
      CALL ErrorStop( global,ERR_NO_PROCMATCH,__LINE__ )

    iRegBeg = 1
    iRegEnd = global%nRegionsProc

    DO iproc=0,global%nProcAlloc-1
      local = 0
      DO iReg=iRegBeg,iRegEnd      ! regions on iproc
        local = local + 1
        IF (global%nProcAlloc>1 .AND. local>=MPI_PATCHOFF) &
          CALL ErrorStop( global,ERR_PATCH_OFFSET,__LINE__ )
        IF (global%nProcAlloc > 1) THEN
          regions(iReg)%procid = iproc
        ELSE
          regions(iReg)%procid = MASTERPROC
        ENDIF
        regions(iReg)%localNumber = local
      ENDDO
      iRegBeg = iRegEnd + 1
      iRegEnd = iRegBeg + global%nRegionsProc - 1
    ENDDO

  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MapRegionsProcessors

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_MapRegionsProcessors.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.6  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.5  2003/02/14 22:32:37  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.2  2002/03/21 18:07:15  jblazek
! Added check of MPI_PATCHOFF (for tags).
!
! Revision 1.1  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.4  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.3  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.2  2001/12/08 00:18:42  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************








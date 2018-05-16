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
! Purpose: write flow solution to file for one region (only mixture for now).
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions            = dimensions and cons. variables of all regions
!        iReg               = region number
!        global%currentTime = physical time
!        global%resInit     = initial residual.
!
! Output: to file.
!
! Notes: solution is stored only for the current grid level; it is also
!        stored for all dummy cells. All regions are written into one file.
!        There is no transfer of data from other processors.
!
!******************************************************************************
!
! $Id: RFLO_WriteSolutionRegion.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_WriteSolutionRegion( iReg,regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_WriteDataFileInt, RFLO_WriteDataFileReal, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: i, j, k, n

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname

  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER :: nDimC, nDimN, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), POINTER     :: cv(:,:), siVel(:), sjVel(:), skVel(:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), cvFile(:,:), sVelFile(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_WriteSolutionRegion',&
  'RFLO_WriteSolutionRegion.F90' )

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(5,1),stat=errorFlag )
  ALLOCATE( rvar(2,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open solution file (only if iReg=1) -----------------------------------------

  IF (iReg == 1) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.sola_', &
                                   global%currentTime
        OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.solb_', &
                                   global%currentTime
        OPEN(IF_SOLUT,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
      rvar(1,1) = global%currentTime
      rvar(2,1) = 1._RFREAL

! - steady flow

    ELSE
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%outDir)//TRIM(global%casename)//'.sola_', &
                                global%currentIter
        OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%outDir)//TRIM(global%casename)//'.solb_', &
                                global%currentIter
        OPEN(IF_SOLUT,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
      rvar(1,1) = 0._RFREAL
      rvar(2,1) = global%resInit
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! - write time and initial residual to file

    CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat,2,1,rvar )

  ENDIF   ! 1st region

! write solution data ---------------------------------------------------------
! get dimensions and pointers

  iLev = regions(iReg)%currLevel
  CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( regions(iReg),iLev,iOff,ijOff )
  ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
  ijkEnd = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
  nDimC  = ijkEnd - ijkBeg + 1

  CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
  nDimN = (regions(iReg)%levels(iLev)%grid%ipc+1)* &
          (regions(iReg)%levels(iLev)%grid%jpc+1)* &
          (regions(iReg)%levels(iLev)%grid%kpc+1)

! allocate memory for data field

  ALLOCATE( cvFile(CV_MIXT_NEQS,nDimC),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  IF (regions(iReg)%mixtInput%moveGrid) THEN
    ALLOCATE( sVelFile(3,nDimN),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ENDIF

! copy solution into data structure

  cv => regions(iReg)%levels(iLev)%mixt%cv

  n = 0
  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        n   = n + 1
        ijk = IndIJK(i,j,k,iOff,ijOff)
        cvFile(1,n) = cv(CV_MIXT_DENS,ijk)
        cvFile(2,n) = cv(CV_MIXT_XMOM,ijk)
        cvFile(3,n) = cv(CV_MIXT_YMOM,ijk)
        cvFile(4,n) = cv(CV_MIXT_ZMOM,ijk)
        cvFile(5,n) = cv(CV_MIXT_ENER,ijk)
      ENDDO
    ENDDO
  ENDDO

  IF (regions(iReg)%mixtInput%moveGrid) THEN
    siVel => regions(iReg)%levels(iLev)%grid%siVel
    sjVel => regions(iReg)%levels(iLev)%grid%sjVel
    skVel => regions(iReg)%levels(iLev)%grid%skVel
    n = 0
    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend
          n   = n + 1
          ijk = IndIJK(i,j,k,iNOff,ijNOff)
          sVelFile(1,n) = siVel(ijk)
          sVelFile(2,n) = sjVel(ijk)
          sVelFile(3,n) = skVel(ijk)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

! write region number and dimensions

  ivar(1,1) = iReg
  ivar(2,1) = regions(iReg)%levels(iLev)%grid%ipc
  ivar(3,1) = regions(iReg)%levels(iLev)%grid%jpc
  ivar(4,1) = regions(iReg)%levels(ilev)%grid%kpc
  ivar(5,1) = regions(iReg)%nDumCells
  CALL RFLO_WriteDataFileInt( global,IF_SOLUT,global%solutFormat,5,1,ivar )

! write data

  CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                               CV_MIXT_NEQS,nDimC,cvFile )

  IF (regions(iReg)%mixtInput%moveGrid) THEN
    CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                 3,nDimN,sVelFile )
  ENDIF

  DEALLOCATE( cvFile,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  IF (regions(iReg)%mixtInput%moveGrid) THEN
    DEALLOCATE( sVelFile,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  ENDIF

! finalize --------------------------------------------------------------------

  IF (iReg == global%nRegions) THEN
    CLOSE(IF_SOLUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_WriteSolutionRegion

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_WriteSolutionRegion.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.6  2002/10/23 18:43:24  jblazek
! Changed temporary pointer arrays into allocatable arrays
! in grid and solution I/O routines.
!
! Revision 1.5  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.4  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.3  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/08/29 21:52:21  jblazek
! Added I/O of grid speeds.
!
! Revision 1.1  2002/06/07 16:40:36  jblazek
! Grid & solution for all regions in one file.
!
!******************************************************************************








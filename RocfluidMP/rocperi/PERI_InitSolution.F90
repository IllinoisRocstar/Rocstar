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
! Purpose: calls routines for initialisation of PERI solutions
!
! Description: initial solution depends on the periodic flow case selected
!
! Input: iReg   = index of current region
!        region = data of current region
!
! Output: solution to start simulation
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_InitSolution.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_InitSolution( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE ModPeriodic, ONLY        : t_peri_input
#ifdef RFLO
  USE ModInterfaces, ONLY      : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE PERI_ModInterfaces, ONLY : PERI_CoCprInitSolution, PERI_CoCnlInitSolution
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER      :: iReg

! ... loop variables
  INTEGER      :: ir, j, k

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER     :: global
  TYPE(t_peri_input), POINTER :: input

  INTEGER      :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend 
  INTEGER      :: iNOff, ijNOff, ijkN
  REAL(RFREAL) :: maxY, minY, maxZ, minZ, sndMax, rcvMax, sndMin, rcvMin
  REAL(RFREAL) :: delta, pidel
  REAL(RFREAL), POINTER :: xyz(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_InitSolution.F90,v $'

  global => regions(iReg)%global
  CALL RegisterFunction( global,'PERI_InitSolution',&
  'PERI_InitSolution.F90' )

! get initial solution depending on kind of periodic flow ---------------------

  IF ((regions(iReg)%periInput%flowKind == PERI_FLOW_CPR) .OR. &
      (regions(iReg)%periInput%flowKind == PERI_FLOW_CHANNEL)) THEN

! - compute channel halfwidth

    maxY = -1000000._RFREAL
    minY =  1000000._RFREAL
    maxZ = -1000000._RFREAL
    minZ =  1000000._RFREAL

#ifdef RFLO
    DO ir = 1,global%nRegions
      IF (regions(ir)%procid==global%myProcid .AND. & ! region active and
          regions(ir)%active==ACTIVE) THEN            ! on my processor

        iLev = regions(ir)%currLevel
        CALL RFLO_GetDimensPhysNodes( regions(ir),iLev,ipnbeg,ipnend, &
                                      jpnbeg,jpnend,kpnbeg,kpnend )
        CALL RFLO_GetNodeOffset( regions(ir),iLev,iNOff,ijNOff )

        xyz  => regions(ir)%levels(iLev)%grid%xyz
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)

        xyz  => regions(ir)%grid%xyz
#endif

        input => regions(ir)%periInput

        IF (input%flowKind==PERI_FLOW_CPR .OR. &
            input%flowKind==PERI_FLOW_CHANNEL) THEN
#ifdef RFLO
          DO j = jpnbeg, jpnend
            ijkN = IndIJK(1 ,j ,1 ,iNOff,ijNOff)
            maxY = MAX( maxY, xyz(YCOORD,ijkN) ) 
            minY = MIN( minY, xyz(YCOORD,ijkN) ) 
          ENDDO

          DO k = kpnbeg, kpnend
            ijkN = IndIJK(1 ,1 ,k ,iNOff,ijNOff)
            maxZ = MAX( maxZ, xyz(ZCOORD,ijkN) ) 
            minZ = MIN( minZ, xyz(ZCOORD,ijkN) ) 
          ENDDO
#endif
#ifdef RFLU
          DO ijkN = 1, regions(ir)%grid%nVert
            maxY = MAX( maxY, xyz(YCOORD,ijkN) ) 
            minY = MIN( minY, xyz(YCOORD,ijkN) ) 
            maxZ = MAX( maxZ, xyz(ZCOORD,ijkN) ) 
            minZ = MIN( minZ, xyz(ZCOORD,ijkN) ) 
          ENDDO
#endif
        ENDIF ! flowKind
#ifdef RFLO
      ENDIF  ! region active and my processor
#endif
    ENDDO   ! ir

#ifndef MPI
    delta = 0.5_RFREAL*(maxY - minY)  ! channel half width
    IF (ABS( global%refLength - delta ) > &
        MAX( global%refLength,1._RFREAL )*PERI_REAL_SMALL) THEN
      CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
           'maxY-minY in the region /= length scale delta' )
    ENDIF
#endif 
#ifdef MPI
    sndMax = maxY
    sndMin = minY
    CALL MPI_ALLREDUCE( sndMax,rcvMax,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                      global%mpiComm, global%mpierr )
    CALL MPI_ALLREDUCE( sndMin,rcvMin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                      global%mpiComm, global%mpierr )
    maxY = rcvMax
    minY = rcvMin

    sndMax = maxZ
    sndMin = minZ
    CALL MPI_ALLREDUCE( sndMax,rcvMax,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                      global%mpiComm, global%mpierr )
    CALL MPI_ALLREDUCE( sndMin,rcvMin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                      global%mpiComm, global%mpierr )
    maxZ = rcvMax
    minZ = rcvMin
#endif

! - save minY and maxY in permanent data and get y and z dimensions
    input => regions(iReg)%periInput
    input%minmax(1) = minY    
    input%minmax(2) = maxY    
    delta = 0.5_RFREAL*(maxY - minY)  ! channel half width
    pidel = maxZ - minZ               ! channel span 

! - check delta
    IF (ABS(global%refLength - delta) > &
            global%refLength*PERI_REAL_SMALL) THEN
      CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
      'length scale delta is not consistent with the grid' )
    ENDIF 

  ENDIF

  IF (regions(iReg)%periInput%flowKind == PERI_FLOW_CPR) THEN

! - check span dimension
    IF (ABS( global%pi*global%refLength - pidel ) > &
        MAX( global%refLength,1._RFREAL )*PERI_REAL_SMALL) THEN
      CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
                     'spanwise dimension is not consistent with the grid' )
    ENDIF 

    CALL PERI_CoCprInitSolution( regions(iReg) )

  ELSEIF (regions(iReg)%periInput%flowKind == PERI_FLOW_CHANNEL) THEN

! - check span dimension
    IF (ABS( global%pi*global%refLength - pidel ) > &
        MAX( global%refLength,1._RFREAL )*PERI_REAL_SMALL) THEN
      CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
                     'spanwise dimension is not consistent with the grid' )
    ENDIF 

    CALL PERI_CoCnlInitSolution( regions(iReg) )

  ELSEIF (regions(iReg)%periInput%flowKind == PERI_FLOW_BOLA) THEN

!    CALL PERI_BolaInitSolution( regions(iReg) )

  ENDIF
  
  global%moduleVar(1) = regions(iReg)%periInput%meanPgrad

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_InitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_InitSolution.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/06/09 01:09:24  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.6  2003/10/20 00:42:06  wasistho
! modified dimension check
!
! Revision 1.5  2003/10/17 20:27:26  wasistho
! fixed the computation of delta
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/04/03 20:58:54  wasistho
! normal extent in one region ifnot MPI
!
! Revision 1.2  2003/04/03 00:35:17  wasistho
! include channel grid check
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************








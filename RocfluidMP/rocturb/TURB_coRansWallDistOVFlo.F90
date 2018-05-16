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
! Purpose: compute distance of each cell centers to the nearest no-slip wall
!
! Description: distance to nearest wall (walldist) using this direct method
!              is computed by collecting coordinates of no-slip wall patches 
!              in all regions and stored in a buffer. The wall distance is 
!              the minimum of the distance between the coordinates in the 
!              buffer and each cell centroid. In this way, global wall distance 
!              is obtained, assuming open view (OV) from cell centers to 
!              the nearest wall.
!
! Input: region = data of current region.
!
! Output: region%levels%turb%lens = turbulence length scale.
!
! Notes: this method is computationally intensive, not suitable for moving
!        grid. More efficient method is provided by other routines.        
!
!******************************************************************************
!
! $Id: TURB_coRansWallDistOVFlo.F90,v 1.7 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CoRansWallDistOV( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_RansWallDistOVPatch, &
                                 TURB_FloExtrapIntCellScal 
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k, n, iPatch, iLev

! ... local variables
  TYPE(t_global), POINTER :: global

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLevBeg,iLevEnd, ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend
  INTEGER :: iCOff,ijCOff,iNOff,ijNOff, iC,ijkN,ijkNi,ijkNj,ijkNk
  INTEGER :: sndDim, rcvDim, ncount, icount, AccWallDim
  INTEGER, POINTER :: wallDim(:)
  REAL(RFREAL) :: distance, xdist, ydist, zdist, locLens, cDES
  REAL(RFREAL) :: adsi, adsj, adsk, dsi(3), dsj(3), dsk(3)
  REAL(RFREAL), POINTER :: xyz(:,:), cofg(:,:), lens(:)

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_CoRansWallDistOV',&
  'TURB_coRansWallDistOVFlo.F90' )

  IF (global%turbWallDim == 0) GOTO 888

! allocation and initialisation -----------------------------------------------

  ALLOCATE( wallDim(global%nRegions) )  ! wall dimension/region, all processors

  IF (.NOT. (global%turbWorkUnused .eqv. .true.)) &
  CALL ErrorStop( global,ERR_TURB_WORKSPACE,__LINE__ ) 
  global%turbWorkDim = NDIR*global%turbWallDim
  ALLOCATE( global%turbWork1D(global%turbWorkDim) )

  global%turbWorkUnused = .FALSE. 
  wallDim(:) = 0
  ncount     = 0
  icount     = 0
  AccWallDim = 0
  global%turbWorkDim = 0

! search for ns patch at finest level (iLev=1) and store wall xyz in
! global%turbWork1D

  DO iReg=1,global%nRegions

    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      DO iPatch=1,regions(iReg)%nPatches
        CALL TURB_RansWallDistOVPatch( regions(iReg), &
                  regions(iReg)%levels(1)%patches(iPatch) )
      ENDDO ! iPatch
      wallDim(iReg) = global%turbWorkDim - AccWallDim
      AccWallDim    = AccWallDim + wallDim(iReg)
    ENDIF

#ifdef MPI
! - distribut wall dimension acquired from each region/processor to all procs.

    CALL MPI_BCAST( wallDim(iReg),1,MPI_INTEGER,regions(iReg)%procId, &
                    global%mpiComm,global%mpierr )
    IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

! - screen output check
!    write(*,*)global%myProcid,iReg,global%turbWorkDim,wallDim(iReg), &
!          NDIR*global%turbWallDim

! - collect wall coordinates (xyz) from all processors to Master processor

    IF (global%myProcid == MASTERPROC) THEN

      IF (wallDim(iReg) > 0) THEN
        IF (regions(iReg)%procid /= global%myProcid) THEN

          CALL MPI_RECV( global%turbWork1D(1+ncount),wallDim(iReg),MPI_RFREAL, &
               regions(iReg)%procId,iReg,global%mpiComm,status,global%mpierr )
          IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

        ENDIF

        ncount = ncount + wallDim(iReg)

! ----- screen output check        
!        write(*,*) 'ncount',ireg,ncount

      ENDIF
    ELSE   
      IF (regions(iReg)%procid == global%myProcid) THEN
        IF (wallDim(iReg) > 0) THEN

! ------- screen output check
!          write(*,*)'send=',ireg

          CALL MPI_SEND( global%turbWork1D(1+icount),wallDim(iReg),MPI_RFREAL, &
                         MASTERPROC,iReg,global%mpiComm,status,global%mpierr )
          IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

          icount = icount + wallDim(iReg)

        ENDIF ! wallDim
      ENDIF   ! region%procId
    ENDIF     ! myProcId
#endif
  ENDDO   ! iReg

! check consistency of total wall dimension

  IF (global%myProcid == MASTERPROC) THEN 
#ifdef MPI
    global%turbWorkDim = ncount
#endif
    IF (global%turbWorkDim /= NDIR*global%turbWallDim) &
    CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
                 'work array dim. is not consistent with wall array dim.')
  ENDIF

! distribute wall coordinates from Master processor to all

#ifdef MPI
  CALL MPI_BCAST( global%turbWork1D,NDIR*global%turbWallDim,MPI_RFREAL, &
                  MASTERPROC,global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

888  CONTINUE

! compute turbulence length scale for RaNS and modified if DES is selected

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      IF (regions(iReg)%mixtInput%moveGrid .eqv. .true.) THEN
        iLevBeg = regions(iReg)%currLevel
        iLevEnd = regions(iReg)%currLevel
      ELSE
        iLevBeg = 1
        iLevEnd = regions(iReg)%nGridLevels
      ENDIF

      cDES = regions(iReg)%turbInput%cDes
 
      DO iLev = iLevBeg, iLevEnd
        CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                                 jpcbeg,jpcend,kpcbeg,kpcend )
        CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
        CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

        xyz  => regions(iReg)%levels(iLev)%grid%xyz
        cofg => regions(iReg)%levels(iLev)%grid%cofg
        lens => regions(iReg)%levels(iLev)%turb%lens 
        lens =  global%refLength

! ----- RaNS wall distance providing no-slip wall exist

        DO k=kpcbeg,kpcend
          DO j=jpcbeg,jpcend
            DO i=ipcbeg,ipcend
              iC       = IndIJK(i,j,k,iCOff,ijCOff) 
              lens(iC) = 1.E+32_RFREAL
              DO n = 1,NDIR*global%turbWallDim,NDIR
                xdist     = cofg(XCOORD,iC)-global%turbWork1D(n)
                ydist     = cofg(YCOORD,iC)-global%turbWork1D(n+1)
                zdist     = cofg(ZCOORD,iC)-global%turbWork1D(n+2)
                distance  = SQRT( xdist**2 + ydist**2 + zdist**2 )
                lens(iC)  = MIN( lens(iC),distance )
              ENDDO ! n

! ----------- screen output check
!              IF(i==30 .and.k==2) write(*,*)ireg,j,lens(iC)
            ENDDO   ! i
          ENDDO     ! j
        ENDDO       ! k 

! ----- DES length scale and flows without no-slip walls

        IF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
            (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA) .OR. &
            (global%turbWallDim == 0)) THEN
          DO k=kpcbeg,kpcend
            DO j=jpcbeg,jpcend
              DO i=ipcbeg,ipcend
                iC       = IndIJK(i  ,j  ,k  ,iCOff,ijCOff) 
                ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
                ijkNi    = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
                ijkNj    = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
                ijkNk    = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
                dsi(1:3) = xyz(XCOORD:ZCOORD,ijkNi)-xyz(XCOORD:ZCOORD,ijkN)
                dsj(1:3) = xyz(XCOORD:ZCOORD,ijkNj)-xyz(XCOORD:ZCOORD,ijkN)
                dsk(1:3) = xyz(XCOORD:ZCOORD,ijkNk)-xyz(XCOORD:ZCOORD,ijkN)
                adsi     = SQRT( dsi(1)*dsi(1)+dsi(2)*dsi(2)+dsi(3)*dsi(3) )
                adsj     = SQRT( dsj(1)*dsj(1)+dsj(2)*dsj(2)+dsj(3)*dsj(3) )
                adsk     = SQRT( dsk(1)*dsk(1)+dsk(2)*dsk(2)+dsk(3)*dsk(3) )
                locLens  = MAX( adsi,adsj,adsk )
                lens(iC) = MIN( lens(iC),cDES*locLens )
              ENDDO   ! i
            ENDDO     ! j
          ENDDO       ! k 
        ENDIF         ! DES or noWall
      ENDDO           ! iLev

! --- extrapolate solution to dummy cells

      CALL TURB_FloExtrapIntCellScal( regions(iReg),lens )

    ENDIF             ! region on this processor and active
  ENDDO               ! iReg

! finalize --------------------------------------------------------------------

  DEALLOCATE( global%turbWork1D, wallDim )
  global%turbWorkUnused = .TRUE.

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CoRansWallDistOV

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_coRansWallDistOVFlo.F90,v $
! Revision 1.7  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.6  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/09 06:36:11  wasistho
! incorporated HDESSA
!
! Revision 1.3  2004/03/13 03:12:16  wasistho
! get rid of flo/flu identifier in TURB_Co.. routines
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:33:31  wasistho
! changed turb nomenclature
!
! Revision 1.8  2004/02/13 03:31:27  wasistho
! fixed bug in case a processor has multi regions
!
! Revision 1.7  2004/02/12 03:46:48  wasistho
! filled in RaNS lengthscale in dummy cells
!
! Revision 1.6  2004/01/21 03:47:10  wasistho
! modify k index in screen output check
!
! Revision 1.5  2003/10/20 00:45:42  wasistho
! initiate lens to reference length
!
! Revision 1.4  2003/10/09 20:48:29  wasistho
! added DES lengthscale coefficient CDES
!
! Revision 1.3  2003/10/07 23:55:30  wasistho
! bug fixed missing nodeOffsets
!
! Revision 1.2  2003/10/07 20:32:16  wasistho
! turbWork2D to turbWork1D
!
!
!******************************************************************************








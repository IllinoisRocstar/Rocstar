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
! $Id: TURB_coRansWallDistOVFlu.F90,v 1.7 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CoRansWallDistOV( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_RansWallDistOVPatch
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iC, n, iPatch

! ... local variables
  TYPE(t_global), POINTER :: global

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iRegGlob, ncount, icount, AccWallDim
  INTEGER, POINTER :: wallDim(:)
  REAL(RFREAL) :: oneThird, distance, xdist, ydist, zdist, locLens, cDES
  REAL(RFREAL), POINTER :: xyz(:,:), cofg(:,:), lens(:)

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_CoRansWallDistOV',&
  'TURB_coRansWallDistOVFlu.F90' )

  IF (global%turbWallDim == 0) GOTO 888

! allocation and initialisation -----------------------------------------------

  ALLOCATE( wallDim(global%nRegions) )  ! wall dimension/region, all processors

  IF (.NOT. (global%turbWorkUnused .eqv. .true.)) &
  CALL ErrorStop( global,ERR_TURB_WORKSPACE,__LINE__ ) 
  global%turbWorkDim = NDIR*global%turbWallDim
  ALLOCATE( global%turbWork1D(global%turbWorkDim) )

  global%turbWorkUnused = .FALSE. 
  oneThird   = 1._RFREAL/3._RFREAL
  wallDim(:) = 0
  ncount     = 0
  icount     = 0
  AccWallDim = 0
  global%turbWorkDim = 0

! search for ns patch and store wall xyz in global%turbWork1D

  DO iReg=1,global%nRegionsLocal
    iRegGlob = regions(iReg)%iRegionGlobal

    DO iPatch=1,regions(iReg)%grid%nPatches
      CALL TURB_RansWallDistOVPatch( regions(iReg), &
                                     regions(iReg)%patches(iPatch) )
    ENDDO ! iPatch
    wallDim(iRegGlob) = global%turbWorkDim - AccWallDim
    AccWallDim        = AccWallDim + wallDim(iRegGlob)

#ifdef MPI
! - distribut wall dimension acquired from each region/processor to all procs.

!c    CALL MPI_BCAST( wallDim(iRegGlob),1,MPI_INTEGER,regions(iReg)%procId, &
!c                    global%mpiComm,global%mpierr )
    IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

! - screen output check
!    write(*,*)global%myProcid,iReg,global%turbWorkDim,wallDim(iReg), &
!          NDIR*global%turbWallDim

! - collect wall coordinates (xyz) from all processors to Master processor

    IF (global%myProcid == MASTERPROC) THEN

      IF (wallDim(iRegGlob) > 0) THEN
!c        IF (regions(iReg)%procid /= global%myProcid) THEN

!c          CALL MPI_RECV( global%turbWork1D(1+ncount),wallDim(iRegGlob), &
!c                         MPI_RFREAL,regions(iReg)%procId,iRegGlob, &
!c                         global%mpiComm,status,global%mpierr )
          IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

!c        ENDIF

        ncount = ncount + wallDim(iRegGlob)

! ----- screen output check        
!        write(*,*) 'ncount',ireg,ncount

      ENDIF
    ELSE   
      IF (wallDim(iRegGlob) > 0) THEN

! ----- screen output check
!       write(*,*)'send=',ireg

        CALL MPI_SEND( global%turbWork1D(1+icount),wallDim(iRegGlob),MPI_RFREAL, &
                       MASTERPROC,iRegGlob,global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

        icount = icount + wallDim(iRegGlob)

      ENDIF ! wallDim
    ENDIF   ! myProcId
#endif
  ENDDO     ! iReg

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

  DO iReg=1,global%nRegionsLocal

    cDES = regions(iReg)%turbInput%cDes

    xyz  => regions(iReg)%grid%xyz
    cofg => regions(iReg)%grid%cofg
    lens => regions(iReg)%turb%lens 
    lens =  global%refLength

! - RaNS wall distance providing no-slip wall exist

    DO iC = 1,regions(iReg)%grid%nCellsTot
      lens(iC) = 1.E+32_RFREAL
      DO n = 1,NDIR*global%turbWallDim,NDIR
        xdist     = cofg(XCOORD,iC)-global%turbWork1D(n)
        ydist     = cofg(YCOORD,iC)-global%turbWork1D(n+1)
        zdist     = cofg(ZCOORD,iC)-global%turbWork1D(n+2)
        distance  = SQRT( xdist**2 + ydist**2 + zdist**2 )
        lens(iC)  = MIN( lens(iC),distance )
      ENDDO ! n

! --- screen output check
!     IF(i==30 .and.k==2) write(*,*)ireg,j,lens(iC)
    ENDDO   ! iC

! - DES length scale and flows without no-slip walls

    IF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
        (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA) .OR. &
        (global%turbWallDim == 0)) THEN
      DO iC = 1,regions(iReg)%grid%nCellsTot
        loclens  = regions(iReg)%grid%vol(iC)**oneThird  ! waiting for 
                                                         ! Andreas` c2e
        lens(iC) = MIN( lens(iC),cDES*locLens )
      ENDDO ! iC
    ENDIF   ! DES or no-nswall

  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  DEALLOCATE( global%turbWork1D, wallDim )
  global%turbWorkUnused = .TRUE.

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CoRansWallDistOV

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_coRansWallDistOVFlu.F90,v $
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
! Revision 1.4  2005/03/09 06:36:15  wasistho
! incorporated HDESSA
!
! Revision 1.3  2005/01/12 01:13:12  wasistho
! removed single quote signs since SUN has trouble with it
!
! Revision 1.2  2004/05/28 01:55:38  wasistho
! commented MPI lines temporarily to compile RFLU with MPI
!
! Revision 1.1  2004/03/25 04:42:57  wasistho
! prepared for RFLU
!
!
!
!******************************************************************************








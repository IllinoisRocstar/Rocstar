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
! Purpose: convert a 2-D grid into ROCFLO`s 3-D format.
!
! Description: none.
!
! Input: none.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TO3D_Main.F90,v 1.4 2008/12/06 08:44:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCFLO_2Dto3D

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TO3D_ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
                                 RFLO_GetNodeOffset, RFLO_WriteGridRegion
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... loop variables
  INTEGER :: i, j, k, ijk, ijkm1

! ... local variables
  CHARACTER(CHRLEN)   :: versionString, headerString
  CHARACTER(3*CHRLEN) :: fname2d
 
  INTEGER :: inOff, ijnOff, ijkBeg, ijkEnd
  INTEGER :: idnbeg, jdnbeg, kdnbeg, idnend, jdnend, kdnend
  INTEGER :: ipnbeg, jpnbeg, kpnbeg, ipnend, jpnend, kpnend
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53
 
  REAL(RFREAL)          :: depth, dz
  REAL(RFREAL), POINTER :: xyz(:,:)
 
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCFLO_2Dto3D',&
  'TO3D_Main.F90' )

! initialize ------------------------------------------------------------------
! global parameters

  global%verbLevel = VERBOSE_HIGH

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = 0._RFREAL    ! no physical time set
  global%currentIter = 0            ! no iteration

  global%inDir  = './'              ! directory path
  global%outDir = './'

  global%nProcAlloc = 1
  global%myProcid   = MASTERPROC    ! default process number (if not MPI)
  global%mpierr     = ERR_NONE
  global%error      = ERR_NONE

  global%pi  = 4._RFREAL*ATAN(1._RFREAL)
  global%rad = global%pi/180._RFREAL

! print header ----------------------------------------------------------------

#ifdef MPI
  CALL MPI_Init( global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  CALL BuildVersionString( versionString )

  headerString = ' '
  versionWidth = LEN_TRIM(versionString)
  margin       = (headerWidth-versionWidth)/2
  headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
  headerString(1:1) = '*'
  headerString(headerWidth:headerWidth) = '*'

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//'  *****************************************************'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *       ROCFLO-MP: 2-D to 3-D Grid Convertor        *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *       ====================================        *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! region data

  global%nRegions = 1

  ALLOCATE( regions(global%nRegions),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  regions(1)%global            => global
  regions(1)%active             = ACTIVE
  regions(1)%procid             = global%myProcid
  regions(1)%nPatches           = 1
  regions(1)%nGridLevels        = 1
  regions(1)%nDumCells          = 0
  regions(1)%mixtInput%moveGrid = .false.  ! fixed grid

  ALLOCATE( regions(1)%levels(regions(1)%nGridLevels),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /=0 ) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! read user input -------------------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Name of file with 2-D grid:'
  READ(STDIN,'(A)') fname2d

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Case name:'
  READ(STDIN,'(A)') global%casename

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Number of cells in the 3rd direction (>1):'
  READ(STDIN,*) regions(1)%levels(1)%grid%kpc
  regions(1)%levels(1)%grid%kpc = MAX(2,regions(1)%levels(1)%grid%kpc)

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Depth of the 3-D grid (>0.):'
  READ(STDIN,*) depth

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Format of 3-D grid (0=ASCII, 1=binary):'
  READ(STDIN,*) global%gridFormat

! read the 2-D grid -----------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading 2-D grid ...'

  OPEN(IF_GRID,FILE=fname2d,status='old',form='formatted',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname2d) )

! dimensions

  READ(IF_GRID,*) regions(1)%levels(1)%grid%ipc, &
                  regions(1)%levels(1)%grid%jpc

  WRITE(STDOUT,'(A,I4,A,I4)') SOLVER_NAME//'   dimensions= ', &
                              regions(1)%levels(1)%grid%ipc, &
                              ' x ',regions(1)%levels(1)%grid%jpc

  CALL RFLO_GetDimensDummyNodes( regions(1),1,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( regions(1),1,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( regions(1),1,inOff,ijnOff )

  ijkBeg = IndIJK(idnbeg,jdnbeg,kdnbeg,inOff,ijnOff)
  ijkEnd = IndIJK(idnend,jdnend,kdnend,inOff,ijnOff)

! allocate memory for coordinates

  ALLOCATE( regions(1)%levels(1)%grid%xyz(3,ijkBeg:ijkEnd),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /=0 ) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! coordinates (store as k=kpnbeg plane at z=depth)

  xyz => regions(1)%levels(1)%grid%xyz

  k = kpnbeg
  DO j=jpnbeg,jpnend
    DO i=ipnbeg,ipnend
      ijk = IndIJK(i,j,k,inOff,ijnOff)
      READ(IF_GRID,*) xyz(XCOORD,ijk), &
                      xyz(YCOORD,ijk)
                      xyz(ZCOORD,ijk) = 0._RFREAL
    ENDDO
  ENDDO

! close file

  CLOSE(IF_GRID,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname2d) )

! translate the x-y-plane -----------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Generating 3-D grid ...'

  dz = depth/REAL(regions(1)%levels(1)%grid%kpc)

  DO k=kpnbeg+1,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk   = IndIJK(i,j,k  ,inOff,ijnOff)
        ijkm1 = IndIJK(i,j,k-1,inOff,ijnOff)
        xyz(XCOORD,ijk) = xyz(XCOORD,ijkm1)
        xyz(YCOORD,ijk) = xyz(YCOORD,ijkm1)
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijkm1) + dz
      ENDDO
    ENDDO
  ENDDO

! store 3-D grid --------------------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Storing 3-D grid ...'

  CALL RFLO_WriteGridRegion( 1,regions )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

END PROGRAM ROCFLO_2Dto3D

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TO3D_Main.F90,v $
! Revision 1.4  2008/12/06 08:44:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:37:55  wasistho
! rflo_modinterfacesto3d to to3d_modinterfaces
!
! Revision 1.1  2004/12/03 02:50:44  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:52:40  wasistho
! lower to upper case
!
! Revision 1.15  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.14  2003/03/20 22:37:45  haselbac
! Renamed ModInterfaces
!
! Revision 1.13  2003/03/20 19:49:34  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.12  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.11  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.10  2002/09/05 17:40:23  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.8  2002/06/14 17:48:36  jblazek
! Added call to MPI_Finalize.
!
! Revision 1.7  2002/06/14 17:46:51  jblazek
! Added version string.
!
! Revision 1.6  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.5  2002/03/18 21:31:55  jblazek
! Made codes compatible with MPI.
!
! Revision 1.4  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.1  2001/12/21 23:56:52  jblazek
! Added utility to convert 2D grids to 3D.
!
!******************************************************************************








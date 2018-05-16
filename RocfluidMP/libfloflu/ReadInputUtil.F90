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
! Purpose: read in a section of a file (until # is encountered), read
!          keywords and store the associated numerical values.
!
! Description: file contains the following subroutines:
!
!  - ReadSection           = section applies to all regions (reads reals)
!  - ReadStringSection     = section applies to all regions (reads strings)
!  - ReadBothSection       = section applies to all regions (reals and strings)
!  - ReadRegionSection     = section applies to a range of regions
!                            (brbeg:brend), and reads reals
!  - ReadBothRegionSection = section applies to a range of regions
!                            (brbeg:brend), and reads both reals and strings
!  - ReadPatchSection      = section applies to a range of patches (prbeg:prend)
!                            within a range of regions (brbeg:brend)
!  - ReadListSection       = section contains a list of values below the keyword
!  - ReadPrefixedListSection = section contains a list of values below the
!                              keyword; each list is preceded by a string
!
! Input: fileID = file number
!        nvals  = number of values to search for and to store
!        keys   = keywords to search for
!        nCols  = no. of columns in the list
!
! Output: vals    = values associated with keywords (reals only)
!         defined = flag if for certain keyword a value was read in
!         brbeg   = begin of region range (values set for these regions)
!         brend   = end of region range
!         prbeg   = begin of patch range (values set for these patches)
!         prend   = end of patch range
!         distrib = single value for a patch (=0) or distribution (>0)
!         fname   = file with distribution for a patch
!         nRows   = no. of rows of the list
!
! Notes:
!
!   In the routines ReadStringSection, ReadBothSection, and
!   ReadBothRegionSection a string does not include the possible
!   comment given to it in the input deck (comments begin with "!").
!
!******************************************************************************
!
! $Id: ReadInputUtil.F90,v 1.5 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadSection( global,fileID,nvals,keys,vals,defined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER      :: fileID, nvals
  CHARACTER(*) :: keys(nvals)
  LOGICAL      :: defined(nvals)
  REAL(RFREAL) :: vals(nvals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc

!******************************************************************************

  CALL RegisterFunction( global,'ReadSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  defined(:) = .false.

  DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT

    DO ival=1,nvals
      IF (.NOT. (defined(ival) .eqv. .true.)) THEN
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          READ(line(nc+1:256),*) vals(ival)
          defined(ival) = .true.
          EXIT
        ENDIF ! line
      ENDIF   ! defined
    ENDDO     ! ival
  ENDDO       ! <empty>

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadSection

! #############################################################################
! #############################################################################

SUBROUTINE ReadStringSection( global,fileID,nvals,keys,vals,defined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: fileID, nvals
  CHARACTER(*)   :: keys(nvals)
  LOGICAL        :: defined(nvals)
  CHARACTER(*)   :: vals(nvals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc, iEnd

!******************************************************************************

  CALL RegisterFunction( global,'ReadStringSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  defined(:) = .false.

  DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT

    DO ival=1,nvals
      IF (.NOT. (defined(ival).eqv..true.)) THEN
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          iEnd = INDEX(line,'!')-1    ! last character before comment begins
          IF (iEnd < 0) iEnd = 256    ! if no comment, retain entire line
          vals(ival)=ADJUSTL(line(nc+1:iEnd))
          defined(ival) = .true.
          EXIT
        ENDIF ! line
      ENDIF   ! defined
    ENDDO     ! ival
  ENDDO       ! <empty>

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadStringSection

! #############################################################################
! #############################################################################

SUBROUTINE ReadBothSection( global,fileID,nvals,nStrVals,keys,strKeys, &
                            vals,strVals,defined,strDefined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER      :: fileID, nvals, nStrVals
  CHARACTER(*) :: keys(nvals), strKeys(nStrVals)
  LOGICAL      :: defined(nvals), strDefined(nStrVals)
  REAL(RFREAL) :: vals(nvals)
  CHARACTER(*) :: strVals(nStrVals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc, iEnd

!******************************************************************************

  CALL RegisterFunction( global,'ReadBothSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  defined(:)    = .false. ! keeps track of values being provided by the user
  strDefined(:) = .false. ! keeps track of string values

o:DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT o

    DO ival=1,nvals
      IF (.NOT. (defined(ival).eqv..true.)) THEN
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          READ(line(nc+1:256),*) vals(ival)
          defined(ival) = .true.
          CYCLE o
        ENDIF ! line
      ENDIF   ! defined
    ENDDO     ! ival

    DO ival=1,nStrVals
      IF (.NOT. (strDefined(ival).eqv..true.)) THEN
        nc = LEN_TRIM(strKeys(ival))
        IF (line(1:nc) == TRIM(strKeys(ival))) THEN ! found matching keyword
          iEnd = INDEX(line,'!')-1    ! last character before comment begins
          IF (iEnd < 0) iEnd = 256    ! if no comment, retain entire line
          strVals(ival)=ADJUSTL(line(nc+1:iEnd))
          strDefined(ival) = .true.
          CYCLE o
        ENDIF ! line
      ENDIF   ! strDefined
    ENDDO     ! ival

  ENDDO o

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadBothSection

! #############################################################################
! #############################################################################

SUBROUTINE ReadRegionSection( global,fileID,nvals,keys,vals, &
                              brbeg,brend,defined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: fileID, nvals, brbeg, brend
  CHARACTER(*)   :: keys(nvals)
  LOGICAL        :: defined(nvals)
  REAL(RFREAL)   :: vals(nvals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc

!******************************************************************************

  CALL RegisterFunction( global,'ReadRegionSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  brbeg = 1               ! region range: input applies to all regions (default)
  brend = global%nRegions

  defined(:) = .false.    ! keeps track of values being provided by the user

  DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT

    IF (line(1:5) == 'BLOCK') THEN
      READ(line(6:256),*) brbeg,brend
      brend = MIN(brend,global%nRegions)
      IF (brbeg <= 0    ) brbeg = 1
      IF (brend <= 0    ) brend = global%nRegions
      IF (brend <  brbeg) brend = brbeg
    ELSE
      DO ival=1,nvals
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          READ(line(nc+1:256),*) vals(ival)
          defined(ival) = .true.
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadRegionSection

! #############################################################################
! #############################################################################

SUBROUTINE ReadBothRegionSection( global,fileID,nvals,nStrVals,keys,strKeys, &
                                  vals,strVals,brbeg,brend,defined,strDefined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER      :: fileID, nvals, nStrVals, brbeg, brend
  CHARACTER(*) :: keys(nvals), strKeys(nStrVals)
  LOGICAL      :: defined(nvals), strDefined(nStrVals)
  REAL(RFREAL) :: vals(nvals)
  CHARACTER(*) :: strVals(nStrVals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc, iEnd

!******************************************************************************

  CALL RegisterFunction( global,'ReadBothRegionSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  brbeg = 1               ! region range: input applies to all regions (default)
  brend = global%nRegions

  defined(:)    = .false. ! keeps track of values being provided by the user
  strDefined(:) = .false. ! keeps track of string values

o:DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT o

    IF (line(1:5) == 'BLOCK') THEN
      READ(line(6:256),*) brbeg,brend
      brend = MIN(brend,global%nRegions)
      IF (brbeg <= 0    ) brbeg = 1
      IF (brend <= 0    ) brend = global%nRegions
      IF (brend <  brbeg) brend = brbeg
      CYCLE o
    ENDIF ! line

    DO ival=1,nvals
      IF (.NOT. (defined(ival).eqv..true.)) THEN
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          READ(line(nc+1:256),*) vals(ival)
          defined(ival) = .true.
          CYCLE o
        ENDIF ! line
      ENDIF   ! defined
    ENDDO     ! ival

    DO ival=1,nStrVals
      IF (.NOT. (strDefined(ival).eqv..true.)) THEN
        nc = LEN_TRIM(strKeys(ival))
        IF (line(1:nc) == TRIM(strKeys(ival))) THEN ! found matching keyword
          iEnd = INDEX(line,'!')-1    ! last character before comment begins
          IF (iEnd < 0) iEnd = 256    ! if no comment, retain entire line
          strVals(ival)=ADJUSTL(line(nc+1:iEnd))
          strDefined(ival) = .true.
          CYCLE o
        ENDIF ! line
      ENDIF   ! strDefined
    ENDDO     ! ival

  ENDDO o

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadBothRegionSection

! #############################################################################
! #############################################################################

#ifdef RFLO
SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals,brbeg,brend, &
                             prbeg,prend,distrib,fname,defined )
#endif

#ifdef RFLU
SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals, &
                             prbeg,prend,distrib,fname,bcName,defined )
#endif


  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  INTEGER      :: brbeg, brend
#endif
  INTEGER      :: fileID, nvals, prbeg, prend, distrib
  CHARACTER(*) :: keys(nvals), fname
#ifdef RFLU
  CHARACTER(*) :: bcName
#endif
  LOGICAL      :: defined(nvals)
  REAL(RFREAL) :: vals(nvals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc

!******************************************************************************

  CALL RegisterFunction( global,'ReadPatchSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

#ifdef RFLO
  brbeg = 1               ! region range: input applies to all regions (default)
  brend = global%nRegions
#endif

  prbeg = 1               ! patch range: input applies to all patches (default)
  prend = 999999          ! can have different # of patches in each region

  distrib    = 0          ! no distribution as a default
  fname      = ''         ! no file name

  IF ( nvals /= 0 ) THEN
    defined(:) = .false.  ! keeps track of values being provided by the user
  END IF ! nvals

#ifdef RFLU
  bcName = 'None'
#endif

  DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT

#ifdef RFLO
    IF (line(1:5) == 'BLOCK') THEN
      READ(line(6:256),*) brbeg,brend
      brend = MIN(brend,global%nRegions)
      IF (brbeg <= 0    ) brbeg = 1
      IF (brend <= 0    ) brend = global%nRegions
      IF (brend <  brbeg) brend = brbeg
    ELSE IF (line(1:5) == 'PATCH') THEN
#endif
#ifdef RFLU
    IF (line(1:5) == 'PATCH') THEN
#endif
      READ(line(6:256),*) prbeg,prend
      IF (prbeg <= 0    ) prbeg = 1
      IF (prend <= 0    ) prend = 999999
      IF (prend <  prbeg) prend = prbeg
    ELSE IF (line(1:7) == 'DISTRIB') THEN
      READ(line(8:256),*) distrib
      distrib = MAX(distrib,0)
      distrib = MIN(distrib,1)
#ifdef RFLO
    ELSE IF (line(1:4) == 'FILE') THEN
      READ(line(5:256),*) fname
#endif
#ifdef RFLU
    ELSE IF (line(1:4) == 'NAME') THEN
      READ(line(5:CHRLEN),*) bcName
      bcName = ADJUSTL(bcName)
    ELSE IF (line(1:4) == 'FILE') THEN
      fname = ADJUSTL(line(5:CHRLEN))
#endif
    ELSE
      DO ival=1,nvals
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          READ(line(nc+1:256),*) vals(ival)
          defined(ival) = .true.
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadPatchSection

! #############################################################################
! #############################################################################

SUBROUTINE ReadListSection( global,fileID,key,nCols,nRows,vals,defined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER      :: fileID, nCols, nRows
  CHARACTER(*) :: key
  LOGICAL      :: defined
  REAL(RFREAL), POINTER :: vals(:,:)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival, n

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc

!******************************************************************************

  CALL RegisterFunction( global,'ReadListSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  defined = .false.    ! initial status
  nRows   = 0

  nc = LEN_TRIM(key)

  DO
    READ(fileID,'(A256)',err=10,end=10) line
    IF (line(1:1) == '#') EXIT

    IF (line(1:nc) == TRIM(key)) THEN
      READ(line(nc+1:256),*,err=10,end=10) nRows
      IF (nRows > 0) THEN
        ALLOCATE( vals(nRows,nCols),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        vals(:,:) = 0.0_RFREAL
        DO ival=1,nRows
          READ(fileID,*,err=10,end=10) (vals(ival,n), n=1,nCols)
        ENDDO
        defined = .true.
      ENDIF
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )

999  CONTINUE

END SUBROUTINE ReadListSection

! #############################################################################
! #############################################################################

SUBROUTINE ReadPrefixedListSection( global,fileID,key,nCols,nRows, &
                                    vals,strVals,defined )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER      :: fileID, nCols, nRows
  CHARACTER(*) :: key
  LOGICAL      :: defined
  REAL(RFREAL), POINTER :: vals(:,:)
  CHARACTER(*), POINTER :: strVals(:)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival, n

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc

!******************************************************************************

  CALL RegisterFunction( global,'ReadListSection',&
  'ReadInputUtil.F90' )

! read lines from file until # or EOF found

  defined = .false.    ! initial status
  nRows   = 0

  nc = LEN_TRIM(key)

  DO
    READ(fileID,'(A256)',err=10,end=10) line
    IF (line(1:1) == '#') EXIT

    IF (line(1:nc) == TRIM(key)) THEN
      READ(line(nc+1:256),*,err=10,end=10) nRows
      IF (nRows > 0) THEN
        ALLOCATE( vals(nRows,nCols),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        vals(:,:) = 0.0_RFREAL
        ALLOCATE( strVals(nRows),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        strVals(:) = ""
        DO ival=1,nRows
          READ(fileID,*,err=10,end=10) strVals(ival),(vals(ival,n), n=1,nCols)
          strVals(ival) = ADJUSTL(strVals(ival))
        ENDDO
        defined = .true.
      ENDIF
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )

999  CONTINUE

END SUBROUTINE ReadPrefixedListSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadInputUtil.F90,v $
! Revision 1.5  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.2  2005/05/03 03:27:25  wasistho
! modified RFLO reading fname
!
! Revision 1.1  2004/12/01 16:50:31  haselbac
! Initial revision after changing case
!
! Revision 1.21  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.20  2004/01/29 22:52:45  haselbac
! Added RFLU support for FILE string
!
! Revision 1.19  2003/12/04 03:23:03  haselbac
! Cosmetic change
!
! Revision 1.18  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.17  2003/03/24 23:25:48  jferry
! moved some routines from libfloflu to rocinteract
!
! Revision 1.16  2003/03/11 22:50:56  jferry
! Added ReadPrefixedListSection subroutine
!
! Revision 1.15  2003/03/11 16:04:19  jferry
! Added ReadBothSection and ReadBothRegionSection subroutines
!
! Revision 1.14  2002/10/08 15:48:35  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.13  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.12  2002/09/10 20:32:48  haselbac
! Temporary workaround for SGI compiler bug
!
! Revision 1.11  2002/09/09 14:03:08  haselbac
! Initialize vals, otherwise get FPE with pgf90
!
! Revision 1.10  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/07/25 14:48:51  haselbac
! Fixed bug for nvals=0
!
! Revision 1.8  2002/06/14 21:17:01  wasistho
! Added time avg statistics
!
! Revision 1.7  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.6  2002/03/26 19:05:53  haselbac
! Added ROCFLU functionality, extended ReadSection to deal with similar names
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.3  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************















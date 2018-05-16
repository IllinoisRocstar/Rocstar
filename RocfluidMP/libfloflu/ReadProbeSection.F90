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
! Purpose: read in user input related to position of a probe.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = location of probe, dump intervall.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadProbeSection.F90,v 1.5 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadProbeSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadListSection, ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival, n

! ... local variables
  CHARACTER(10) :: keys(3)  
  LOGICAL :: defined(3)
  INTEGER :: errorFlag, nCols, nRows  
  REAL(RFREAL) :: valsDump(3)
  REAL(RFREAL), POINTER :: valsLoc(:,:)

!******************************************************************************

  CALL RegisterFunction( global,'ReadProbeSection',&
  'ReadProbeSection.F90' )

! do not read probes a second time

! TEMPORARY - Will be fixed properly later, when we will have routines to 
!             create and destroy probes. Error trapping apparently only needed
!             because reading this again will cause allocation to be executed
!             again. This is a problem for rflumap.
!  IF ( global%nProbes > 0 ) THEN 
!    CALL ErrorStop(global,ERR_PROBE_SPECIFIED,__LINE__)
!  END IF ! global%nProbes
  
  IF ( global%nProbes == 0 ) THEN 
! END TEMPORARY  

! specify keywords and search for them

#ifdef RFLO
  defined(:) = .false.

  keys(1) = 'NUMBER'
  nCols   = 4
  
  CALL ReadListSection( global,IF_INPUT,keys(1),nCols,nRows,valsLoc,defined(1) )

  IF (defined(1).eqv..true.) THEN
    global%nProbes = nRows
    ALLOCATE( global%probePos(nRows,nCols),stat=errorFlag )
    errorFlag = global%error
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - support entering 0 for the block and the x, y, z coordinates
!   instead of the block number and indeces.

    ALLOCATE( global%probeXYZ(nRows,nCols),stat=errorFlag )
    errorFlag = global%error
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    DO ival=1,global%nProbes
      IF (valsLoc(ival,1) /= 0.) THEN
        DO n=1,nCols
          global%probePos(ival,n) = INT(ABS(valsLoc(ival,n))+0.5_RFREAL)
        ENDDO
      ELSE

! ----- they have entered coordinates.  Assign probePos in writeProbe.
        DO n=1,nCols
          global%probePos(ival,n) = 0
          global%probeXYZ(ival,n) = valsLoc(ival,n)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  IF (defined(1).eqv..true.) THEN
    DEALLOCATE( valsLoc,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  ENDIF

! get dump interval

  defined(:) = .false.

  keys(1) = 'WRITIME'
  keys(2) = 'WRIITER'
  keys(3) = 'OPENCLOSE'

  CALL ReadSection( global,IF_INPUT,3,keys,valsDump,defined )

  IF (defined(1).eqv..true.) global%probeSaveTime = ABS(valsDump(1))
  IF (defined(2).eqv..true.) THEN
    global%probeSaveIter = INT(ABS(valsDump(2))+0.5_RFREAL)
    global%probeSaveIter = MAX(1,global%probeSaveIter)
  ENDIF
  IF (defined(3).eqv..true.) THEN
    IF (valsDump(3) < 0.5_RFREAL) THEN
      global%probeOpenClose = .false.
    ELSE
      global%probeOpenClose = .true.
    ENDIF
  ENDIF
  IF ((.NOT.(defined(1).eqv..true.)).AND. &
      (.NOT.(defined(2).eqv..true.)).AND. &
      (.NOT.(defined(3).eqv..true.))) THEN
    BACKSPACE(IF_INPUT, IOSTAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) &
      CALL ErrorStop( global,ERR_FILE_READ,__LINE__, &
                     'while backspacing after reading probe section' )
  ENDIF  ! not defined
#endif


#ifdef RFLU
  defined = .FALSE.

  keys(1) = 'NUMBER'
  nCols   = 3

  CALL ReadListSection(global,IF_INPUT,keys(1),nCols,nRows,valsLoc,defined(1))

  IF ( defined(1) .EQV. .TRUE. ) THEN
    global%nProbes = nRows

    ALLOCATE(global%probePos(nRows,PROBE_REGION:PROBE_CELL),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)    
    END IF ! global%error
    
    DO ival = 1,global%nProbes
      DO n = PROBE_REGION,PROBE_CELL
        global%probePos(ival,n) = CRAZY_VALUE_INT
      END DO ! n
    END DO ! ival    
    
    ALLOCATE(global%probeXYZ(nRows,nCols),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
    END IF ! global%error

    DO ival = 1,global%nProbes
      DO n = 1,nCols
        global%probeXYZ(ival,n) = valsLoc(ival,n)
      END DO ! n
    END DO ! ival    
  END IF ! defined  

  IF ( defined(1) .EQV. .TRUE. ) THEN
    DEALLOCATE(valsLoc,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
    END IF ! global%error
  END IF ! defined

! get dump interval

  defined(:) = .FALSE.

  keys(1) = 'WRITIME'
  keys(2) = 'WRIITER'
  keys(3) = 'OPENCLOSE'

  CALL ReadSection(global,IF_INPUT,3,keys,valsDump,defined)

  IF ( defined(1) .EQV. .TRUE. ) THEN 
    global%probeSaveTime = ABS(valsDump(1))
  END IF ! defined
  IF ( defined(2) .EQV. .TRUE. ) THEN
    global%probeSaveIter = INT(ABS(valsDump(2))+0.5_RFREAL)
    global%probeSaveIter = MAX(1,global%probeSaveIter)
  END IF ! defined
  IF ( defined(3) .EQV. .TRUE. ) THEN
    IF ( valsDump(3) < 0.5_RFREAL ) THEN
      global%probeOpenClose = .FALSE.
    ELSE
      global%probeOpenClose = .TRUE.
    END IF ! valsDump
  END IF ! defined
#endif

! TEMPORARY - See comment above
  END IF ! global%nProbes
! END TEMPORARY

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadProbeSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadProbeSection.F90,v $
! Revision 1.5  2008/12/06 08:44:10  mtcampbe
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
! Revision 1.2  2006/03/25 02:17:57  wasistho
! added safety when certain params not exist
!
! Revision 1.1  2004/12/01 16:50:44  haselbac
! Initial revision after changing case
!
! Revision 1.17  2004/11/11 14:49:56  haselbac
! Commented out error check for probes, bypass for now
!
! Revision 1.16  2004/08/09 22:14:42  fnajjar
! Changed apostrophe in comment line since SUN compiler breaks
!
! Revision 1.15  2004/07/21 21:11:42  wasistho
! allow probes input by coordinates
!
! Revision 1.14.2.1  2004/07/02 04:09:27  rfiedler
! Allows Rocflo probes to be specified by coordinates.  Use 0 for the block ID.
!
! Revision 1.14  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.11  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.10  2003/04/07 14:18:40  haselbac
! Added new options for RFLU
!
! Revision 1.9  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.8  2002/10/08 15:48:35  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.7  2002/10/05 18:37:11  haselbac
! Added allocation of probeXyz
!
! Revision 1.6  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/03/26 19:07:20  haselbac
! Added ROCFLU functionality
!
! Revision 1.4  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
!******************************************************************************








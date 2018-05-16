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
! Purpose: read in user input for ROCVAV.
!
! Description: none.
!
! Input: user input file.
!
! Output: None
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ReadInputFile.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadInputFile( global ) 

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE ModParameters
  USE RVAV_ModParameters                          
  USE RVAV_ModGlobal
  USE RVAV_ModInterfaces, ONLY : RVAV_ReadSectionStream1, &
                                 RVAV_ReadSectionStream2, &
                                 RVAV_ReadComparisonsSection
  IMPLICIT NONE

! ... parameter variables
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(CHRLEN+4) :: fname
  CHARACTER(256)      :: line

  INTEGER :: errorFlag

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ReadInputFile',&
  'RVAV_ReadInputFile.F90' )

! open file

  WRITE(fname,'(A,I5.5,A,1PE11.5)') TRIM(global%casename)//'.inp_RVAV'
  
  OPEN(IF_RVAV_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global, ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  DO
    READ(IF_RVAV_INPUT,'(A256)',err=10,end=999) line
        
    SELECT CASE(TRIM(line))
    CASE ('# RVAV_STREAM1')
      IF ( global%verbLevel/=VERBOSE_NONE ) &
        WRITE(STDOUT,'(/,A)') '     RVAV_ReadSectionStream1 ...'
      CALL RVAV_ReadSectionStream1( global ) ! reads everything specified about stream1
                                             ! in the input file RVAV_Input.inp

    CASE ('# RVAV_STREAM2')
      IF ( global%verbLevel/=VERBOSE_NONE ) &
        WRITE(STDOUT,'(/,A)') '     RVAV_ReadStream2Section ...'
      CALL RVAV_ReadSectionStream2( global ) ! reads everything specified about stream2
                                             ! in the input file RVAV_Input.inp

    CASE ('# RVAV_COMPARISONS')
      IF ( global%verbLevel/=VERBOSE_NONE ) &
        WRITE(STDOUT,'(/,A)') '     RVAV_ReadComparisonSection ...'
      CALL RVAV_ReadComparisonsSection( global ) ! reads everything about the comparison
                                                 ! to be made with the available data

! if there are more details in the input file put in the appropriate
! subroutine calls here 

    END SELECT
  ENDDO
  
  IF ( global%verbLevel/=VERBOSE_NONE ) &
    WRITE(STDOUT,'(/,A)') 'Closing IF_RVAV_INPUT ...'

  CLOSE(IF_RVAV_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

  GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RVAV_ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ReadInputFile.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:29  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.4  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.2  2002/08/15 19:48:06  jblazek
! Implemented grid deformation capability.
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************








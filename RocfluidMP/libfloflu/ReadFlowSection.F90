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
! Purpose: read in user input related to flow model.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = flow model, moving grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadFlowSection.F90,v 1.4 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadFlowSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
#ifdef RFLO  
  USE ModInterfaces, ONLY : ReadRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : ReadSection
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(10)     :: keys(2)
 
#ifdef RFLO 
  INTEGER :: brbeg, brend
#endif 
#ifdef RFLU
  INTEGER :: iReg
#endif 
 
  LOGICAL :: defined(2)
 
  REAL(RFREAL) :: vals(2)

!******************************************************************************

  CALL RegisterFunction( regions(1)%global,'ReadFlowSection',&
  'ReadFlowSection.F90' )

! specify keywords and search for them

  keys(1) = 'MODEL'
  keys(2) = 'MOVEGRID'

#ifdef RFLO
  CALL ReadRegionSection( regions(1)%global,IF_INPUT,2,keys,vals, &
                          brbeg,brend,defined )

  IF (defined(1).eqv..true.) THEN
                     regions(brbeg:brend)%mixtInput%flowModel = FLOW_EULER
    IF (vals(1)>0.9) regions(brbeg:brend)%mixtInput%flowModel = FLOW_NAVST
  ENDIF
  IF (defined(2).eqv..true.) THEN
    IF (vals(2) < 0.1) THEN
      regions(brbeg:brend)%mixtInput%moveGrid = .false.
    ELSE
      regions(brbeg:brend)%mixtInput%moveGrid = .true.
    ENDIF
  ENDIF 
#endif

#ifdef RFLU
  CALL ReadSection( regions(1)%global,IF_INPUT,2,keys,vals,defined ) 
  
  IF ( defined(1) .eqv..true.) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%flowModel = NINT(vals(1))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(2) .eqv..true.) THEN
    IF ( NINT(vals(2)) == 0 ) THEN
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
        regions(iReg)%mixtInput%moveGrid = .FALSE.
      END DO ! iReg
    ELSEIF ( NINT(vals(2)) == 1 ) THEN      
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)
        regions(iReg)%mixtInput%moveGrid = .TRUE.
      END DO ! iReg
    END IF ! NINT
  END IF ! defined 
#endif 

! finalize

  CALL DeregisterFunction( regions(1)%global )

END SUBROUTINE ReadFlowSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadFlowSection.F90,v $
! Revision 1.4  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.1  2004/12/01 16:50:12  haselbac
! Initial revision after changing case
!
! Revision 1.16  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.12  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.11  2003/03/15 16:22:59  haselbac
! Added KIND qualifyers
!
! Revision 1.10  2002/10/16 21:09:43  haselbac
! Fixed bug in RFLU code segment
!
! Revision 1.9  2002/09/09 14:01:16  haselbac
! mixtInput now under regions
!
! Revision 1.8  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.7  2002/03/26 19:01:46  haselbac
! Added ROCFLU functionality
!
! Revision 1.6  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.4  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
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








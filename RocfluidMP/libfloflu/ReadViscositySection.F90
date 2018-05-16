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
! Purpose: read in user input related to viscosity model.
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
! $Id: ReadViscositySection.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadViscositySection( regions )

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
  INTEGER, PARAMETER :: NKEYS = 4
  CHARACTER(10) :: keys(NKEYS)
 
#ifdef RFLO 
  INTEGER :: brbeg, brend
#endif 
#ifdef RFLU
  INTEGER :: iReg
#endif 
 
  LOGICAL :: defined(NKEYS)
 
  REAL(RFREAL) :: vals(NKEYS)

!******************************************************************************

  CALL RegisterFunction( regions(1)%global,'ReadViscositySection',&
  'ReadViscositySection.F90' )

! specify keywords and search for them

  keys(1) = 'MODEL'
  keys(2) = 'VISCOSITY'
  keys(3) = 'REFTEMP'
  keys(4) = 'SUTHCOEF'

#ifdef RFLO
  CALL ReadRegionSection( regions(1)%global,IF_INPUT,NKEYS,keys,vals, &
                          brbeg,brend,defined )

  IF (defined(1).eqv..true.) &
    regions(brbeg:brend)%mixtInput%viscModel = NINT(vals(1))

  IF (defined(2).eqv..true.) THEN
    regions(brbeg:brend)%mixtInput%refVisc   = ABS(vals(2))
  ELSE
    regions(brbeg:brend)%mixtInput%refVisc   = -1._RFREAL
  ENDIF
  
  IF (defined(3).eqv..true.) &
    regions(brbeg:brend)%mixtInput%refTemp   = ABS(vals(3))
  
  IF (defined(4).eqv..true.) &
    regions(brbeg:brend)%mixtInput%suthCoef  = ABS(vals(4))    
#endif

#ifdef RFLU
  CALL ReadSection( regions(1)%global,IF_INPUT,NKEYS,keys,vals,defined ) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%viscModel = NINT(vals(1))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%refVisc = ABS(vals(2))
    END DO ! iReg
  ELSE 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%refVisc = REAL(CRAZY_VALUE_INT,RFREAL)
    END DO ! iReg    
  END IF ! defined 
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%refTemp = ABS(vals(3))
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%suthCoef = ABS(vals(4))
    END DO ! iReg
  END IF ! defined 
#endif 

! finalize

  CALL DeregisterFunction( regions(1)%global )

END SUBROUTINE ReadViscositySection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadViscositySection.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
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
! Revision 1.1  2004/12/01 16:50:59  haselbac
! Initial revision after changing case
!
! Revision 1.8  2003/12/04 03:23:04  haselbac
! Added parameter, setting of refVisc if undefined, fixed bug
!
! Revision 1.7  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/20 22:16:16  jblazek
! Corrected bug in viscosity model input.
!
! Revision 1.3  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/04/10 23:28:27  fnajjar
! Initial import for viscosity models
!
!******************************************************************************








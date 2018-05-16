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
! Purpose: Read in user input within PERI section (done on all processors).
!
! Description: none.
!
! Input: regions = user input file of all regions.
!
! Output: regions = PERI input parameters.
!
! Notes: Mother routine = ReadInputFile.
!
!******************************************************************************
!
! $Id: PERI_ReadPeriSection.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_ReadPeriSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO 
  USE ModInterfaces, ONLY : ReadRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : ReadSection
#endif
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global


  INTEGER :: nVals
  INTEGER :: brbeg, brend
  INTEGER, PARAMETER :: NVALS_MAX = 6

  LOGICAL           :: defined(NVALS_MAX)
  REAL(RFREAL)      :: vals(NVALS_MAX)
  CHARACTER(15)     :: keys(NVALS_MAX)

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_ReadPeriSection.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'PERI_ReadPeriSection',&
  'PERI_ReadPeriSection.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX
  keys( 1) = 'FLOWKIND'
  keys( 2) = 'ISPLIT'
  keys( 3) = 'JSPLIT'
  keys( 4) = 'KSPLIT'
  keys( 5) = 'CPREPSILON'
  keys( 6) = 'PGRADTYPE'

#ifdef RFLO

  CALL ReadRegionSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                          brbeg,brend,defined(1:nVals) )
#endif
#ifdef RFLU  
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) ) 

  brbeg = LBOUND(regions,1)   ! temporary for now before zonal modeling apply
  brend = UBOUND(regions,1)
#endif

  IF (defined(1)) THEN
    regions(brbeg:brend)%periInput%flowKind = INT(vals(1)+0.5_RFREAL)
  ENDIF
  IF (defined(2)) THEN
    regions(brbeg:brend)%periInput%split(ICOORD) = INT(vals(2)+0.5_RFREAL)
  ENDIF
  IF (defined(3)) THEN
    regions(brbeg:brend)%periInput%split(JCOORD) = INT(vals(3)+0.5_RFREAL)
  ENDIF
  IF (defined(4)) THEN
    regions(brbeg:brend)%periInput%split(KCOORD) = INT(vals(4)+0.5_RFREAL)
  ENDIF
  IF (defined(5)) THEN
    regions(brbeg:brend)%periInput%cprEpsilon = ABS(vals(5))
  ENDIF
  IF (defined(6)) THEN
    regions(brbeg:brend)%periInput%pgradType = INT(vals(6)+0.5_RFREAL)
  ENDIF

! finalize -------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_ReadPeriSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_ReadPeriSection.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.3  2003/09/18 01:57:04  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.2  2003/04/02 01:48:22  wasistho
! minimize CPR user input
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************








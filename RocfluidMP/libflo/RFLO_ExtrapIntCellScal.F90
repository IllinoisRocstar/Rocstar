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
! Purpose: Extrapolate scalar variable defined at cells to outer dummy cells.
!
! Description: The values are only available in the interior cells.
!              Here the values are extrapolated to first up to nDum cells.
!
! Input: region = data of current region
!        fVec   = scalar variable to be extrapolated
!
! Output: fVec = values of fVec at the remaining dummy cells are defined.
!
!******************************************************************************
!
! $Id: RFLO_ExtrapIntCellScal.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExtrapIntCellScal( region,fVec )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy,RFLO_GetCellOffset
  USE ModTurbulence
  USE ModError
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region)        :: region
  REAL(RFREAL), POINTER :: fVec(:)

! ... loop variables
  INTEGER :: i, j, k, l, m, ijkC, ijkC0, ijkC1, ijkC2, ijkC3
  INTEGER :: indxC, indxC0, indxC1, indxC2, indxC3

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER :: iLev,iCOff,ijCOff,nDum,nadd,sign(2),indx(2)
  INTEGER :: iDum

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLO_ExtrapIntCellScal.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'RFLO_ExtrapIntCellScal',&
  'RFLO_ExtrapIntCellScal.F90' )

! get indices ---------------------------------------------------------------

  sign(1)= 1
  sign(2)=-1
  nDum   = region%nDumCells
  nadd   = -nDum
  iLev   = region%currLevel

  CALL RFLO_GetDimensDummy( region,ilev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

  ibeg = idcbeg
  iend = idcend
  jbeg = jdcbeg
  jend = jdcend
  kbeg = kdcbeg
  kend = kdcend

! extrapolation to dummy layers, first in I direction

  indx(1)=ibeg
  indx(2)=iend
  DO m=1,2
    DO iDum=2,nDum 
      indxC3 = indx(m)-sign(m)*(nadd-iDum+1)
      indxC2 = indx(m)-sign(m)*nadd+sign(m)
      indxC1 = indx(m)-sign(m)*nadd 
      indxC0 = indx(m)-sign(m)*nadd-sign(m)
      indxC  = indx(m)-sign(m)*(nadd+iDum)
      DO k=kbeg,kend
        DO j=jbeg,jend
          ijkC  = IndIJK(indxC   ,j  ,k  ,iCOff,ijCOff)
          ijkC0 = ijkC + indxC0-indxC
          ijkC1 = ijkC + indxC1-indxC
          ijkC2 = ijkC + indxC2-indxC
          ijkC3 = ijkC + indxC3-indxC

          fVec(ijkC0) = 2._RFREAL*fVec(ijkC1)-fVec(ijkC2) 
          fVec(ijkC)  = fVec(ijkC0)+fVec(ijkC1)-fVec(ijkC3) 
        ENDDO   ! j
      ENDDO     ! k
    ENDDO       ! iDum
  ENDDO         ! m

! then in J direction

  indx(1)=jbeg
  indx(2)=jend
  DO m=1,2
    DO iDum=2,nDum 
      indxC3 = indx(m)-sign(m)*(nadd-iDum+1)
      indxC2 = indx(m)-sign(m)*nadd+sign(m)
      indxC1 = indx(m)-sign(m)*nadd 
      indxC0 = indx(m)-sign(m)*nadd-sign(m)
      indxC  = indx(m)-sign(m)*(nadd+iDum)
      DO k=kbeg,kend
        DO i=ibeg,iend
          ijkC  = IndIJK(i  ,indxC   ,k  ,iCOff,ijCOff)
          ijkC0 = ijkC + (indxC0-indxC)*iCOff
          ijkC1 = ijkC + (indxC1-indxC)*iCOff
          ijkC2 = ijkC + (indxC2-indxC)*iCOff
          ijkC3 = ijkC + (indxC3-indxC)*iCOff

          fVec(ijkC0) = 2._RFREAL*fVec(ijkC1)-fVec(ijkC2) 
          fVec(ijkC)  = fVec(ijkC0)+fVec(ijkC1)-fVec(ijkC3) 
        ENDDO   ! i
      ENDDO     ! k
    ENDDO       ! iDum
  ENDDO         ! m

! then in K direction

  indx(1)=kbeg
  indx(2)=kend
  DO m=1,2
    DO iDum=2,nDum 
      indxC3 = indx(m)-sign(m)*(nadd-iDum+1)
      indxC2 = indx(m)-sign(m)*nadd+sign(m)
      indxC1 = indx(m)-sign(m)*nadd 
      indxC0 = indx(m)-sign(m)*nadd-sign(m)
      indxC  = indx(m)-sign(m)*(nadd+iDum)
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkC  = IndIJK(i  ,j  ,indxC   ,iCOff,ijCOff)
          ijkC0 = ijkC + (indxC0-indxC)*ijCOff
          ijkC1 = ijkC + (indxC1-indxC)*ijCOff
          ijkC2 = ijkC + (indxC2-indxC)*ijCOff
          ijkC3 = ijkC + (indxC3-indxC)*ijCOff

          fVec(ijkC0) = 2._RFREAL*fVec(ijkC1)-fVec(ijkC2) 
          fVec(ijkC)  = fVec(ijkC0)+fVec(ijkC1)-fVec(ijkC3) 
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! iDum
  ENDDO         ! m

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExtrapIntCellScal

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExtrapIntCellScal.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.1  2004/09/30 17:02:39  wasistho
! added RFLO_extrapIntCellVec/Scal routines
!
!
!
!******************************************************************************








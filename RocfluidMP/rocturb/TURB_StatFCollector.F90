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
! Purpose: Collect variables of interest to be time averaged for verification
!
! Description: As a generic routine, number and kind of collected variables
!              are not specified. The cell values are obtained by averaging 
!              the face values, including dummy cells. Data is accumulated
!              in I, J and K direction, respectively, through three subsquent
!              calls.
!
! Input: region = data of current region
!        ijk    = averaging direction
!        iBegSv = begin index of averaged field sv
!        iEndSv = end index of averaged field sv
!        colVar = array of collected variables
!
! Output: turb%st(iBegSv:iEndSv,:) at cell centers including dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_StatFCollector.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE TURB_StatFCollector( region,ijk,iBegSt,iEndSt,colVar )
#endif
#ifdef RFLU
SUBROUTINE TURB_StatFCollector( region,ijk,iBegSt,iEndSt,colVar,colBVar )
#endif

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: ijk, iBegSt, iEndSt
  REAL(RFREAL), POINTER :: colVar(:,:)
#ifdef RFLU
  REAL(RFREAL), POINTER :: colBVar(:,:)  ! collected variables at boundaries
#endif

! ... loop variables
  INTEGER :: i, j, k, l, m

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ijkC, ijkN
  REAL(RFREAL) :: factor
  REAL(RFREAL), POINTER :: st(:,:)

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev,iCOff,ijCOff,iNOff,ijNOff,ijkNP,iadd
#endif
#ifdef RFLU
  INTEGER :: ict, icl, iPatch, nFacesPerCell
  INTEGER, POINTER :: c2f(:,:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_StatFCollector.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'Turb_StatFCollector',&
  'TURB_StatFCollector.F90' )

! check some input arguments --------------------------------------------------

  IF (iEndSt > region%turbInput%nSt) THEN
    CALL ErrorStop( global,ERR_TURB_STATSINPUT,__LINE__, &
        'index of collected vars larger than nSt (allocated nmbr stats vars)' )
  ENDIF
  IF (iBegSt > iEndSt) THEN
    CALL ErrorStop( global,ERR_TURB_STATSINPUT,__LINE__, &
        'begin index of collected vars larger than end index' )
  ENDIF

#ifdef RFLO
! get parameters --------------------------------------------------------

  iLev  = region%currLevel

! get dimensions and pointers

  CALL RFLO_GetDimensDummy( region,ilev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  st => region%levels(ilev)%turb%st

! interpolate from face to cell

  factor = 1._RFREAL/6._RFREAL
  IF (ijk == DIRI) THEN
    st = 0._RFREAL
    iadd = 1
  ELSEIF (ijk == DIRJ) THEN
    iadd = iNOff
  ELSEIF (ijk == DIRK) THEN
    iadd = ijNOff
  ENDIF

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend

        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        ijkNP = ijkN + iadd

        DO l=iBegSt,iEndSt   
          m = l - iBegSt + 1  
          st(l,ijkC) = st(l,ijkC)+factor*(colVar(m,ijkN)+colVar(m,ijkNP))
        ENDDO

      ENDDO ! i
    ENDDO   ! j
  ENDDO     ! k
#endif
#ifdef RFLU
! get dimensions and pointers ------------------------------------------

  st => region%turb%st
  st =  0._RFREAL

  DO ijkC=1,region%grid%nCells

! - region (here denoted as 'global') to local (type) mapping

    ict = region%grid%cellGlob2Loc(1,ijkC) ! cell type
    icl = region%grid%cellGlob2Loc(2,ijkC) ! local (type) cell index
    SELECT CASE ( ict ) 
      CASE ( CELL_TYPE_TET ) 
        nFacesPerCell = 4
        c2f => region%grid%tet2f              
      CASE ( CELL_TYPE_HEX ) 
        nFacesPerCell = 6
        c2f => region%grid%hex2f              
      CASE ( CELL_TYPE_PRI ) 
        nFacesPerCell = 5    
        c2f => region%grid%pri2f                         
      CASE ( CELL_TYPE_PYR ) 
        nFacesPerCell = 5
        c2f => region%grid%pyr2f                
      CASE DEFAULT  
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ict

    DO l=iBegSt,iEndSt   
      m = l - iBegSt + 1  
      DO k=1,nFacesPerCell
        iPatch = c2f(1,k,icl)    ! patch #, ipatch=0 denotes interior faces
        ijkN   = c2f(2,k,icl)    ! global face index

        IF ( iPatch == 0 ) THEN  ! Interior face
          st(l,ijkC) = st(l,ijkC)+colVar(m,ijkN)
        ELSE ! Boundary face
          st(l,ijkC) = st(l,ijkC)+colBVar(m,ijkN)
        ENDIF
      ENDDO
      st(l,ijkC) = st(l,ijkC)/nFacesPerCell
    ENDDO   ! l
  ENDDO     ! ijkC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_StatFCollector

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_StatFCollector.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/10/22 23:19:34  wasistho
! fixed devision by 6 faces
!
! Revision 1.3  2004/05/28 02:01:48  wasistho
! update unstructured grid LES
!
! Revision 1.2  2004/03/23 03:35:00  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2003/05/24 02:30:28  wasistho
! turbulence statistics expanded
!
!
!******************************************************************************








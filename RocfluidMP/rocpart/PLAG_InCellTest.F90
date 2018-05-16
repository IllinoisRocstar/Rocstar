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
! Purpose: Perform in cell test
!
! Description: none.
!
! Input: region = data of current region
!        posPlag = particle position vector
!        indexSearch = indices for cell search
!        ijkNR,ijkNRI,ijkNRJ,ijkNRK = cell indices for all faces.
!        
! Output: indexNew = new cell index
!          cellLocate = logical variable set to TRUE if test successful
!
! Notes: reverse sign of face vectors to point inward for lbound=1,3,5
!        to check inside cell. 
!        use scaled values for face normals to perform test.
!
!******************************************************************************
!
! $Id: PLAG_InCellTest.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_inCellTest(region, posPlag, indexSearch, &
                           ijkNR,ijkNRI,ijkNRJ,ijkNRK,   &
                           indexNew,cellLocate)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError  
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER                   , INTENT(IN)    :: ijkNR,ijkNRI,ijkNRJ,ijkNRK
  INTEGER,      DIMENSION(4), INTENT(IN)    :: indexSearch
  INTEGER,      DIMENSION(4), INTENT(OUT)   :: indexNew
  
  LOGICAL,                    INTENT(OUT)   :: cellLocate  
  
  REAL(RFREAL), DIMENSION(3), INTENT(IN)    :: posPlag  

! ... loop variables
  INTEGER :: lbound
  
! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nBound

  REAL(RFREAL), PARAMETER    :: epsDegenTol = -1.0E-10_RFREAL
  REAL(RFREAL)               :: dpFace
  REAL(RFREAL), DIMENSION(3) :: diffPos, faceCentroid, sFace
  
  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: pSNormal
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pFc
  
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InCellTest.F90,v $ $Revision: 1.3 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_inCellTest',&
  'PLAG_InCellTest.F90' )

! Get dimensions --------------------------------------------------------------

  iLev = region%currLevel
    
  nBound = 6
  cellLocate = .FALSE.

! Set pointers ----------------------------------------------------------------

  pFc => region%levels(iLev)%plag%fc

! Loop over all cell faces ----------------------------------------------------
 
  DO lbound = 1, nBound
  
    SELECT CASE (lbound)

! - i-face check --------------------------------------------------------------
   
      CASE(1)
        sface(1:3)        = -region%levels(iLev)%plag%si(XCOORD:ZCOORD,ijkNR)
        faceCentroid(1:3) = pFc(XCOORD:ZCOORD,ICOORD,ijkNR)

      CASE(2)
        sface(1:3)        =  region%levels(iLev)%plag%si(XCOORD:ZCOORD,ijkNRI)
        faceCentroid(1:3) = pFc(XCOORD:ZCOORD,ICOORD,ijkNRI)

! - j-face check -------------------------------------------------------------- 
   
      CASE(3)
        sface(1:3)        = -region%levels(iLev)%plag%sj(XCOORD:ZCOORD,ijkNR)
        faceCentroid(1:3) = pFc(XCOORD:ZCOORD,JCOORD,ijkNR)

      CASE(4)
        sface(1:3)        =  region%levels(iLev)%plag%sj(XCOORD:ZCOORD,ijkNRJ)
        faceCentroid(1:3) = pFc(XCOORD:ZCOORD,JCOORD,ijkNRJ)

! - k-face check -------------------------------------------------------------  
  
      CASE(5)
        sface(1:3)        = -region%levels(iLev)%plag%sk(XCOORD:ZCOORD,ijkNR)
        faceCentroid(1:3) = pFc(XCOORD:ZCOORD,KCOORD,ijkNR)

      CASE(6)
        sface(1:3)        =  region%levels(iLev)%plag%sk(XCOORD:ZCOORD,ijkNRK)
        faceCentroid(1:3) = pFc(XCOORD:ZCOORD,KCOORD,ijkNRK)

    END SELECT ! lbound

! - Compute position vector difference ----------------------------------------
!   and perform dot product with face vectors ---------------------------------
!   need to check (r_p - r_fc). n_fc > 0 --------------------------------------
  
    diffPos(1:3) = posPlag(1:3)-faceCentroid(1:3)
    dpFace = DOT_PRODUCT( sFace,diffPos )
  
! - exit immediately if particle is outside the cell --------------------------

    IF ( dpFace < epsDegenTol ) GOTO 999

  ENDDO ! lbound 
  
! Particle passed all face tests ----------------------------------------------
  
  cellLocate = .TRUE.
  indexNew(1:4) = indexSearch(1:4)
  
! finalize --------------------------------------------------------------------

999  CONTINUE
  CALL DeregisterFunction( global )
 
END SUBROUTINE PLAG_inCellTest

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InCellTest.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:35  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/03/23 15:54:21  fnajjar
! Defined epsDegenTol for skewed grids allowing tolerance in particle search algorithm
!
! Revision 1.8  2003/11/03 21:21:51  fnajjar
! Changed definition of face vectors pointing to PLAG datastructure
!
! Revision 1.7  2003/04/25 21:09:04  jferry
! various streamlining, including removing scaling of normals by area
!
! Revision 1.6  2003/04/18 19:21:30  fnajjar
! Redefined dpFace to be a true dot product using FORTRAN90 intrinisic
!
! Revision 1.5  2003/04/17 01:31:34  fnajjar
! Scaled normals by area
!
! Revision 1.4  2003/04/17 00:11:23  fnajjar
! Included Proper INTENT for calling sequence
!
! Revision 1.3  2003/04/16 22:33:42  fnajjar
! Bug fix for selecting appropriate normals of face vectors
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








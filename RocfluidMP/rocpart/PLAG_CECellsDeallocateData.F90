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
! Purpose: deallocate memory for variables associated with buffer datastructure 
!          of the corner and edge cells for all active regions.
!
! Description: none.
!
! Input: regions   = all regions,
!        iReg      = region number.
!
! Output: region%level%cornerCells(:)%bufferExchPlag = corner cell buffers
!         region%level%edgeCells(:)%bufferExchPlag   =   edge cell buffers
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsDeallocateData.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_cECellsDeallocateData( regions, iReg )

  USE ModDataTypes 
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag 
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE ModMPI

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: iCorner, iEdge, iLev, ijk 

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: errorFlag, nGridLevels

  TYPE(t_region),      POINTER :: pRegion
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff, pEdgeCellsXBuff
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsDeallocateData.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global, 'PLAG_CECellsDeAllocateData',&
  'PLAG_CECellsDeallocateData.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE     ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,  &
    'Deallocating Corner-Edge Cells Data Buffers for PLAG...'
  END IF ! global%verbLevel

! Set pointer -----------------------------------------------------------------
  
  pRegion => regions(iReg)
 
! Get dimensions --------------------------------------------------------------

  nGridLevels  = pRegion%nGridLevels 
       
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,nGridLevels

! - Set pointers --------------------------------------------------------------
  
    pLevel => pRegion%levels(iLev)
    pPlag  => regions(iReg)%levels(iLev)%plag

! - Corner cells --------------------------------------------------------------
 
    DO iCorner=1,8

! -- Bypass for noninteracting regions ----------------------------------------

        IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 1999

! -- Bypass for degenerate corner cells ---------------------------------------

      IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 1999

! -- Set pointer --------------------------------------------------------------

      DO ijk=1,UBOUND(pLevel%cornerCells(iCorner)%cells,1)
        pCornCellsXBuff => pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag
      
! -- Deallocate buffer data ---------------------------------------------------

        DEALLOCATE( pCornCellsXBuff%aiv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%aiv' ) 
        END IF ! global%error 

        DEALLOCATE( pCornCellsXBuff%arv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%arv' )
        END IF ! global%error

        DEALLOCATE( pCornCellsXBuff%cv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%cv' ) 
        END IF ! global%error

        DEALLOCATE( pCornCellsXBuff%dv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%dv' ) 
        END IF ! global%error

        DEALLOCATE( pCornCellsXBuff%tv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%tv' ) 
        END IF ! global%error

        DEALLOCATE( pCornCellsXBuff%aivOld,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%aivOld' ) 
        END IF ! global%error 

        DEALLOCATE( pCornCellsXBuff%arvOld,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%arvOld' ) 
        END IF ! global%error

        DEALLOCATE( pCornCellsXBuff%cvOld,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%cvOld' ) 
        END IF ! global%error  

        DEALLOCATE( pCornCellsXBuff%rhs,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%rhs' ) 
        END IF ! global%error  

        DEALLOCATE( pCornCellsXBuff%rhsSum,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pCornCellsXBuff%rhsSum' ) 
        END IF ! global%error 
      ENDDO ! ijk
      
1999  CONTINUE

    ENDDO ! iCorner

! - Edge cells ----------------------------------------------------------------

    DO iEdge=1,12

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 2999

! -- Bypass for degenerate edge cells -----------------------------------------

      IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! -- Set pointer --------------------------------------------------------------

      DO ijk=1,UBOUND(pLevel%edgeCells(iedge)%cells,1)
        pEdgeCellsXBuff   => pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag

! -- Deallocate buffer data ---------------------------------------------------

        DEALLOCATE( pEdgeCellsXBuff%aiv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%aiv' ) 
        END IF ! global%error 

        DEALLOCATE( pEdgeCellsXBuff%arv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%arv' )
        END IF ! global%error

        DEALLOCATE( pEdgeCellsXBuff%cv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%cv' ) 
        END IF ! global%error

        DEALLOCATE( pEdgeCellsXBuff%dv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%dv' ) 
        END IF ! global%error

        DEALLOCATE( pEdgeCellsXBuff%tv,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%tv' ) 
        END IF ! global%error

        DEALLOCATE( pEdgeCellsXBuff%aivOld,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%aivOld' ) 
        END IF ! global%error 

        DEALLOCATE( pEdgeCellsXBuff%arvOld,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%arvOld' ) 
        END IF ! global%error

        DEALLOCATE( pEdgeCellsXBuff%cvOld,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%cvOld' ) 
        END IF ! global%error  

        DEALLOCATE( pEdgeCellsXBuff%rhs,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%rhs' ) 
        END IF ! global%error  

        DEALLOCATE( pEdgeCellsXBuff%rhsSum,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__,'pEdgeCellsXBuff%rhsSum' ) 
        END IF ! global%error 
      ENDDO ! ijk

2999  CONTINUE

    ENDDO ! iEdge
    
  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsDeallocateData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsDeallocateData.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:08  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/11/29 19:24:17  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.1  2003/11/12 21:37:59  fnajjar
! Initial import of Corner-Edge cells Infrastructure
!
!******************************************************************************








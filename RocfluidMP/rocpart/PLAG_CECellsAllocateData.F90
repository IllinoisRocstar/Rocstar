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
! Purpose: allocate memory for variables associated with buffer datastructure 
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
! Notes: an improved kernel is to determine the maximum buffer sizes
!        of all the regions and use that value to allocate the arrays.
!
!******************************************************************************
!
! $Id: PLAG_CECellsAllocateData.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_cECellsAllocateData( regions, iReg )

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

  INTEGER :: errorFlag, nAiv, nArv, nCont, nCv, nDv, nTv
  INTEGER :: nBuffI, nBuffR, nBuffSizeMax, nGridLevels
  INTEGER :: nBuffSizeCorn, nBuffSizeEdge
         
  TYPE(t_region),      POINTER :: pRegion
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff, pEdgeCellsXBuff
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsAllocateData.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global, 'PLAG_CECellsAllocateData',&
  'PLAG_CECellsAllocateData.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE     ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,  &
    'Allocating Corner-Edge Cells Data Buffers for PLAG...'
  END IF ! global%verbLevel

! Set pointer -----------------------------------------------------------------
  
  pRegion => regions(iReg)
 
! Get dimensions --------------------------------------------------------------

  nCont        = pRegion%plagInput%nCont
  nBuffSizeMax = pRegion%plagInput%nPclsBuffCECellsMax
  nGridLevels  = pRegion%nGridLevels 
  nBuffSizeCorn = nBuffSizeMax
  nBuffSizeEdge = nBuffSizeMax
       
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,nGridLevels

! - Set pointers --------------------------------------------------------------
  
    pLevel => pRegion%levels(iLev)
    pPlag  => regions(iReg)%levels(iLev)%plag
 
! - Get dimensions ------------------------------------------------------------

    nAiv = pPlag%nAiv
    nArv = pPlag%nArv
            
    nCv  = pPlag%nCv
    nDv  = pPlag%nDv
    nTv  = pPlag%nTv  
    
    nBuffI = 2*nAiv
    nBuffR = 2*nArv +4*nCv +nDv +nTv

! - Corner cells --------------------------------------------------------------

! - Initialize buffer size for all corner cells -------------------------------
 
    DO iCorner=1,8

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 1999

! -- Set pointer --------------------------------------------------------------

      DO ijk=1,UBOUND(pLevel%cornerCells(iCorner)%cells,1)
        pCornCellsXBuff => pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag    
        pCornCellsXBuff%nBuffSize    = 0
        pCornCellsXBuff%nBuffSizeDes = 0
      ENDDO ! ijk
            
1999  CONTINUE
    ENDDO ! iCorner

! - Allocate buffer arrays only for non-degenerate corner cells ---------------
 
    DO iCorner=1,8

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 2999

! -- Bypass for degenerate corner cells ---------------------------------------

      IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! -- Set pointer --------------------------------------------------------------

      DO ijk=1,UBOUND(pLevel%cornerCells(iCorner)%cells,1)
        pCornCellsXBuff => pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag

        pCornCellsXBuff%nBuffSizeTot = nBuffSizeCorn

! -- Allocate buffer data -----------------------------------------------------

        ALLOCATE( pCornCellsXBuff%aiv(nAiv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%aiv' ) 
        END IF ! global%error 

        ALLOCATE( pCornCellsXBuff%arv(nArv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%arv' )
        END IF ! global%error

        ALLOCATE( pCornCellsXBuff%cv(nCv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%cv' ) 
        END IF ! global%error

        ALLOCATE( pCornCellsXBuff%dv(nDv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%dv' ) 
        END IF ! global%error

        ALLOCATE( pCornCellsXBuff%tv(nTv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%tv' ) 
        END IF ! global%error

        ALLOCATE( pCornCellsXBuff%aivOld(nAiv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%aivOld' ) 
        END IF ! global%error 

        ALLOCATE( pCornCellsXBuff%arvOld(nArv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%arvOld' ) 
        END IF ! global%error

        ALLOCATE( pCornCellsXBuff%cvOld(nCv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%cvOld' ) 
        END IF ! global%error  

        ALLOCATE( pCornCellsXBuff%rhs(nCv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%rhs' ) 
        END IF ! global%error  

        ALLOCATE( pCornCellsXBuff%rhsSum(nCv,nBuffSizeCorn),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pCornCellsXBuff%rhsSum' ) 
        END IF ! global%error 
          
! -- Initialize data --------------------------------------------------------

        pCornCellsXBuff%aiv = 0
        pCornCellsXBuff%arv = 0.0_RFREAL
        pCornCellsXBuff%cv  = 0.0_RFREAL
        pCornCellsXBuff%dv  = 0.0_RFREAL
        pCornCellsXBuff%tv  = 0.0_RFREAL

        pCornCellsXBuff%aivOld = 0
        pCornCellsXBuff%arvOld = 0.0_RFREAL
        pCornCellsXBuff%cvOld  = 0.0_RFREAL

        pCornCellsXBuff%rhs    = 0.0_RFREAL
        pCornCellsXBuff%rhsSum = 0.0_RFREAL
      ENDDO ! ijk
            
2999  CONTINUE

    ENDDO ! iCorner

! - Edge cells ----------------------------------------------------------------

! - Initialize buffer size for all corner cells -------------------------------
 
    DO iEdge=1,12

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 3999

! -- Set pointer --------------------------------------------------------------

      DO ijk=1,UBOUND(pLevel%edgeCells(iedge)%cells,1)
        pEdgeCellsXBuff => pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag
        pEdgeCellsXBuff%nBuffSize    = 0
        pEdgeCellsXBuff%nBuffSizeDes = 0
      ENDDO ! ijk 
      
3999  CONTINUE
    ENDDO ! iEdge

! - Initialize buffer size for all corner cells -------------------------------
 
    DO iEdge=1,12

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 4999

! -- Bypass for degenerate edge cells -----------------------------------------

      IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 4999
      
! -- Set pointer --------------------------------------------------------------

      DO ijk=1,UBOUND(pLevel%edgeCells(iedge)%cells,1)
        pEdgeCellsXBuff => pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag

!      PRINT*, 'ASSOCIATED(cells)= ',ASSOCIATED(pLevel%edgeCells(iEdge)%cells)
!      PRINT*, 'ASSOCIATED(pEdgeCellsXBuff)= ',ASSOCIATED(pEdgeCellsXBuff)

        pEdgeCellsXBuff%nBuffSizeTot = nBuffSizeEdge

! -- Allocate buffer data -----------------------------------------------------

        ALLOCATE( pEdgeCellsXBuff%aiv(nAiv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%aiv' ) 
        END IF ! global%error 

        ALLOCATE( pEdgeCellsXBuff%arv(nArv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%arv' )
        END IF ! global%error

        ALLOCATE( pEdgeCellsXBuff%cv(nCv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%cv' ) 
        END IF ! global%error

        ALLOCATE( pEdgeCellsXBuff%dv(nDv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%dv' ) 
        END IF ! global%error

        ALLOCATE( pEdgeCellsXBuff%tv(nTv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%tv' ) 
        END IF ! global%error

        ALLOCATE( pEdgeCellsXBuff%aivOld(nAiv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%aivOld' ) 
        END IF ! global%error 

        ALLOCATE( pEdgeCellsXBuff%arvOld(nArv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%arvOld' ) 
        END IF ! global%error

        ALLOCATE( pEdgeCellsXBuff%cvOld(nCv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%cvOld' ) 
        END IF ! global%error  

        ALLOCATE( pEdgeCellsXBuff%rhs(nCv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%rhs' ) 
        END IF ! global%error  

        ALLOCATE( pEdgeCellsXBuff%rhsSum(nCv,nBuffSizeEdge),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pEdgeCellsXBuff%rhsSum' ) 
        END IF ! global%error 
          
! -- Initialize data --------------------------------------------------------

        pEdgeCellsXBuff%aiv = 0
        pEdgeCellsXBuff%arv = 0.0_RFREAL
        pEdgeCellsXBuff%cv  = 0.0_RFREAL
        pEdgeCellsXBuff%dv  = 0.0_RFREAL
        pEdgeCellsXBuff%tv  = 0.0_RFREAL

        pEdgeCellsXBuff%aivOld = 0
        pEdgeCellsXBuff%arvOld = 0.0_RFREAL
        pEdgeCellsXBuff%cvOld  = 0.0_RFREAL

        pEdgeCellsXBuff%rhs    = 0.0_RFREAL
        pEdgeCellsXBuff%rhsSum = 0.0_RFREAL
      ENDDO ! ijk 
      
4999  CONTINUE

    ENDDO ! iEdge    
  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsAllocateData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsAllocateData.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:05  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/11/29 19:23:11  fnajjar
! Added bypass statement for dengerate cells and initialized buffer sizes
!
! Revision 1.3  2004/03/10 23:10:48  fnajjar
! Allocating corner-edge cells data with nPclsBuffCECellsMax
!
! Revision 1.2  2004/01/28 21:46:28  fnajjar
! Defined nBuffSizeCorn & nBuffSizeEdge for memory allocation of buffers
!
! Revision 1.1  2003/11/12 21:37:59  fnajjar
! Initial import of Corner-Edge cells Infrastructure
!
!******************************************************************************








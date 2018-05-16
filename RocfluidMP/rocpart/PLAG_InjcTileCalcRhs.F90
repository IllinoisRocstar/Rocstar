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
! Purpose: computes the RHS for the multiphase injection algorithm.
!
! Description: none.
!
! Input: region    = current region
!
! Output: regions(iReg)%levels%patch%tile%rhs = updated tile rhs values
!                                               of current region.
!
! Notes: Use negative values of rhs for consistent RK Updating step.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileCalcRhs.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileCalcRhs( region )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset,   RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, iCont, i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, distrib, i2dVals,  ibeg, iCOff, idir, iend,     &
             ijCOff, ijNOff, ijkC, ijkD, ijkN, iLev, inode, iNOff,   &
             iTile, jbeg, jdir, jend, jnode, kbeg, kdir, kend,       &
             knode, lbound, n1, n2, nCont, nOff, nPatches   
  INTEGER, POINTER, DIMENSION(:) :: pCvTileMass

  REAL(RFREAL)  :: area, heatCapSum, injcVelRatio, massFluxSum, mRate, &
                   rhoMixtbCond, tBurn, tileTemp, tileVelNrm
  REAL(RFREAL), POINTER, DIMENSION(:)   :: injcMassFluxRatio, specHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pRhs, sFace, vals
  
  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global),    POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileCalcRhs.F90,v $ $Revision: 1.4 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_InjcTileCalcRhs',&
  'PLAG_InjcTileCalcRhs.F90' )

! Get dimensions --------------------------------------------------------------

  iLev     = region%currLevel
  nPatches = region%nPatches

  nCont        = region%plagInput%nCont
  injcVelRatio = region%plagInput%injcVelRatio
  
  injcMassFluxRatio => region%plagInput%injcMassFluxRatio
  specHeat          => region%plagInput%spht
  
! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    pPatch  => region%levels(iLev)%patches(iPatch)

    bcType = pPatch%bcType

! - Select injection boundary condition ---------------------------------------

    IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN

! -- Get dimensions and pointers ----------------------------------------------

      lbound = pPatch%lbound

      CALL RFLO_GetPatchIndices( region,pPatch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( pPatch,idir,jdir,kdir )
      CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
      CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

      nOff        = ABS(pPatch%l1end-pPatch%l1beg) + 1
      distrib     = pPatch%mixt%distrib
  
      vals => pPatch%mixt%vals
      
! -- Set appropriate extent of cells and shift upper cell by one value --------

      inode = 0
      jnode = 0
      knode = 0
      IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
        inode =  -idir
        jnode =  -jdir
        knode =  -kdir
      ENDIF ! lbound
      
! -- Get the appropriate face vector ------------------------------------------

      IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%plag%si
      IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%plag%sj
      IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%plag%sk

! -- Loop over patch ----------------------------------------------------------

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkC  = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            ijkD  = IndIJK(i-idir,j-jdir,k-kdir,iCOff,ijCOff)
            area  = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                         sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                         sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg + 1
              n2 = k - kbeg + 1
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg + 1
              n2 = i - ibeg + 1
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg + 1
              n2 = j - jbeg + 1
            ENDIF ! lbound
            
            i2dVals = distrib * IndIJ(n1-1,n2-1,nOff)
            mRate   = vals(BCDAT_INJECT_MFRATE,i2dVals)
            tBurn   = vals(BCDAT_INJECT_TEMP  ,i2dVals)
 
! --- Compute mixture density on boundary -------------------------------------
! --- based on an average of interior & dummy cell ----------------------------

            rhoMixtbCond = 0.5_RFREAL *                                     &
                          (region%levels(iLev)%mixt%cv(CV_MIXT_DENS,ijkC) + &
                           region%levels(iLev)%mixt%cv(CV_MIXT_DENS,ijkD) )

! --- Update tile Rhs datastructure -------------------------------------------
                           
            iTile  = IndIJ(n1,n2,nOff)
            pTilePlag => pPatch%tilePlag
 
            pRhs        => pTilePlag%rhs
            pCvTileMass => pTilePlag%cvTileMass

            massFluxSum = SUM ( mRate * injcMassFluxRatio(:) )
            heatCapSum  = DOT_PRODUCT( specHeat, mRate * &
                                       injcMassFluxRatio(:) )

            tileTemp    = tBurn
            tileVelNrm  = injcVelRatio * mRate / rhoMixtbCond
  
            DO iCont = 1, nCont
              pRhs(pCvTileMass(iCont),iTile) = -area *mRate * &
                                                injcMassFluxRatio(iCont) 
            END DO ! iCont            
            
            pRhs(CV_TILE_MOMNRM,iTile) = -area *massFluxSum *tileVelNrm
            
            pRhs(CV_TILE_ENER  ,iTile) = -area *                       &
                                  ( 0.5_RFREAL*massFluxSum*tileVelNrm**2 &
                                  +            heatCapSum *tileTemp      )
           
          END DO ! i
        END DO ! j
      END DO ! k
                                                        
    ENDIF  ! bcType
 
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InjcTileCalcRhs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileCalcRhs.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:05  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/12/01 20:57:43  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/11/03 21:21:51  fnajjar
! Changed definition of face vectors pointing to PLAG datastructure
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/05/07 15:12:10  fnajjar
! Bug fix of i2dVals for correct boundary extent
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************








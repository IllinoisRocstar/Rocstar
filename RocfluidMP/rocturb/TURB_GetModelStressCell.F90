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
! Purpose: Get subgrid or Reynolds stresses at cell centers.
!
! Description: The cell values are obtained by averaging the face values.
!              The result are stored in turb%sv at cells, including dummmies.
!
! Input: region = data of current region
!
! Output: turb%sv(E11:E33,:) at cell centers incl. dummies.
!
! Notes: Dummy stress obtained by extrapolation as dummy mut not computed.
!
!******************************************************************************
!
! $Id: TURB_GetModelStressCell.F90,v 1.8 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_GetModelStressCell( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l, iC, iN, iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  REAL(RFREAL), POINTER :: mueT(:,:),sIij(:,:), sJij(:,:), sKij(:,:), sv(:,:)
#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev,iCOff,ijCOff,iNOff,ijNOff,ijkC,ijkCi,ijkN,ijkNI,ijkNJ,ijkNK
  REAL(RFREAL) :: one6th
#endif
#ifdef RFLU
  INTEGER :: iC0, iC1, ict, icl, ifg, ifgBeg, nFacesPerCell
  INTEGER, POINTER      :: c2f(:,:,:)
  REAL(RFREAL), POINTER :: bMueT(:,:), bsIij(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_GetModelStressCell.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'Turb_GetModelStressCell',&
  'TURB_GetModelStressCell.F90' )

! get parameters -------------------------------------------------------------

#ifdef RFLO
  one6th = 1._RFREAL/6._RFREAL
  iLev   = region%currLevel

! get dimensions and pointers

  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetDimensDummy( region,ilev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  mueT => region%levels(ilev)%turb%mueT
  sIij => region%levels(ilev)%turb%mISij
  sJij => region%levels(ilev)%turb%mJSij
  sKij => region%levels(ilev)%turb%mKSij
  sv   => region%levels(ilev)%turb%sv

! perform averaging of stress tensor

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend

        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        ijkNI = ijkN + 1
        ijkNJ = ijkN + iNOff
        ijkNK = ijkN + ijNOff

        DO l = E11, E33
          sv(l,ijkC) = &
          mueT(DIRI,ijkN)*sIij(l,ijkN)+mueT(DIRI,ijkNI)*sIij(l,ijkNI)+ &
          mueT(DIRJ,ijkN)*sJij(l,ijkN)+mueT(DIRJ,ijkNJ)*sIij(l,ijkNJ)+ &
          mueT(DIRK,ijkN)*sKij(l,ijkN)+mueT(DIRK,ijkNK)*sIij(l,ijkNK)
          sv(l,ijkC) = one6th*sv(l,ijkC)
        ENDDO ! l
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif

#ifdef RFLU
  mueT  => region%turb%mueT
  bMuet => region%turb%bMueT
  sIij  => region%turb%mISij
  bsIij => region%turb%bmISij
  sv    => region%turb%sv
  sv(:,:) = 0._RFREAL

! compute stress tensor at faces then average it from faces to cells 
! waiting for Andreas` routine for proper averaging

  DO iN = 1,region%grid%nFaces
    iC0 = region%grid%f2c(1,iN)
    iC1 = region%grid%f2c(2,iN)
    DO l = E11, E33
      sv(l,iC0) = sv(l,iC0) + mueT(DIRI,iN)*sIij(l,iN)
      sv(l,iC1) = sv(l,iC1) + mueT(DIRI,iN)*sIij(l,iN)
    ENDDO  ! l
  ENDDO    ! iN

  DO iPatch = 1,region%grid%nPatches
! TEMPORARY : removing usage of bf2bg from everywhere
!    ifgBeg = region%patches(iPatch)%bf2bg(BF2BG_BEG)
    DO iN = 1,region%patches(iPatch)%nBFaces
      iC0 = region%patches(iPatch)%bf2c(iN)
      ifg = iN + ifgBeg-1
      DO l = E11, E33
        sv(l,iC0) = sv(l,iC0) + bMueT(DIRI,ifg)*bsIij(l,ifg)
      ENDDO  ! l
    ENDDO    ! iN
  ENDDO      ! iPatch

  DO iC = 1,region%grid%nCells
    ict = region%grid%cellGlob2Loc(1,iC) ! cell type
    icl = region%grid%cellGlob2Loc(2,iC) ! local (type) cell index
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

    DO l = E11, E33
      sv(l,iC) = sv(l,iC)/nFacesPerCell
    ENDDO  ! l
  ENDDO    ! iC

#endif

! linear extrapolate to the six sides of region

#ifdef RFLO
  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend

      ijkN  = IndIJK(ipcbeg   ,j    ,k     ,iNOff,ijNOff)
      DO i=idcbeg,ipcbeg-1
        ijkC  = IndIJK(i            ,j     ,k     ,iCOff,ijCOff)
        ijkCi = IndIJK(2*ipcbeg-i-1 ,j     ,k     ,iCOff,ijCOff)

        DO l = E11, E33
          sv(l,ijkC) = 2._RFREAL*mueT(DIRI,ijkN)*sIij(l,ijkN) - &
                                 sv(l,ijkCi)
        ENDDO       
      ENDDO  
      ijkN  = IndIJK(ipcend+1  ,j    ,k     ,iNOff,ijNOff)
      DO i=ipcend+1,idcend
        ijkC  = IndIJK(i            ,j     ,k     ,iCOff,ijCOff)
        ijkCi = IndIJK(2*ipcend-i+1 ,j     ,k     ,iCOff,ijCOff)

        DO l = E11, E33
          sv(l,ijkC) = 2._RFREAL*mueT(DIRI,ijkN)*sIij(l,ijkN) - &
                                 sv(l,ijkCi)
        ENDDO       
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
  DO k=kdcbeg,kdcend
    DO i=idcbeg,idcend

      ijkN  = IndIJK(i    ,jpcbeg   ,k    ,iNOff,ijNOff)
      DO j=jdcbeg,jpcbeg-1
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkCi = IndIJK(i     ,2*jpcbeg-j-1 ,k     ,iCOff,ijCOff)

        DO l = E11, E33
          sv(l,ijkC) = 2._RFREAL*mueT(DIRJ,ijkN)*sJij(l,ijkN) - &
                                 sv(l,ijkCi)
        ENDDO       
      ENDDO  
      ijkN  = IndIJK(i    ,jpcend+1 ,k    ,iNOff,ijNOff)
      DO j=jpcend+1,jdcend
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkCi = IndIJK(i     ,2*jpcend-j+1 ,k     ,iCOff,ijCOff)

        DO l = E11, E33
          sv(l,ijkC) = 2._RFREAL*mueT(DIRJ,ijkN)*sJij(l,ijkN) - &
                                 sv(l,ijkCi)
        ENDDO       
      ENDDO   ! j
    ENDDO     ! i
  ENDDO       ! k
  DO j=jdcbeg,jdcend
    DO i=idcbeg,idcend

      ijkN  = IndIJK(i    ,j     ,kpcbeg ,iNOff,ijNOff)
      DO k=kdcbeg,kpcbeg-1
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkCi = IndIJK(i     ,j     ,2*kpcbeg-k-1 ,iCOff,ijCOff)

        DO l = E11, E33
          sv(l,ijkC) = 2._RFREAL*mueT(DIRK,ijkN)*sKij(l,ijkN) - &
                                 sv(l,ijkCi)
        ENDDO       
      ENDDO  
      ijkN  = IndIJK(i    ,j    ,kpcend+1 ,iNOff,ijNOff)
      DO k=kpcend+1,kdcend
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkCi = IndIJK(i     ,j     ,2*kpcend-k+1 ,iCOff,ijCOff)

        DO l = E11, E33
          sv(l,ijkC) = 2._RFREAL*mueT(DIRK,ijkN)*sKij(l,ijkN) - &
                                 sv(l,ijkCi)
        ENDDO       
      ENDDO   ! k
    ENDDO     ! i
  ENDDO       ! j
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_GetModelStressCell

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_GetModelStressCell.F90,v $
! Revision 1.8  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/08/19 15:40:55  mparmar
! Removed bf2bg
!
! Revision 1.5  2005/12/29 19:53:05  wasistho
! modified face indexing in f2c averaging
!
! Revision 1.4  2005/01/12 01:12:53  wasistho
! removed single quote signs since SUN has trouble with it
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/19 02:46:48  wasistho
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








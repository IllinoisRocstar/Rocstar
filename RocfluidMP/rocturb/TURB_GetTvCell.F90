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
! Purpose: Get turbulent variables mu_t, kapp_t and cdyn at cell centers.
!
! Description: The cell values are obtained by averaging the face values.
!              Mu_t and kappa_t are defined at all cells, including dummmies.
!              Cdyn, stored in turb%dv, is available only at interior cells.
!
! Input: region  = data of current region
!
! Output: mixt%tv(TV_MIXT_MUET:TV_MIXT_TCOT,:) at cell centers and
!         turb%dv(DV_TURB_CDYN,:) at interior cells.
!
! Notes: Mu_t and kapp_t at dummy cells are obtained by extrapolation.
!        If needed they can be computed directly in routine LesCalcEddyVis.
!
!******************************************************************************
!
! $Id: TURB_GetTvCell.F90,v 1.14 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_GetTvCell( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE TURB_ModInterfaces, ONLY : TURB_FloExtrapolCellVec
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iC, ijkC, iN, ijkN, iPatch, ift, ifl

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: indCp
  REAL(RFREAL), POINTER :: mueT(:,:),tv(:,:), gv(:,:), tdv(:,:)
  REAL(RFREAL)          :: rPrt,cpPrt

#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev,iCOff,ijCOff,iNOff,ijNOff,ijkNI,ijkNJ,ijkNK
  REAL(RFREAL) :: one6th
#endif
#ifdef RFLU
  INTEGER :: ifg, ifgBeg, iC0, iC1, ict, icl, nFacesPerCell
  INTEGER, POINTER      :: c2f(:,:,:)
  REAL(RFREAL), POINTER :: bMueT(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_GetTvCell.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'Turb_GetTvCell',&
  'TURB_GetTvCell.F90' )

! get parameters ----------------------------------------------------

#ifdef RFLO
  iLev   = region%currLevel
  one6th = 1._RFREAL/6._RFREAL
  rPrt   = 1._RFREAL/region%levels(iLev)%mixt%prTurb
  indCp  = region%levels(iLev)%mixt%indCp

! get dimensions and pointers

  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetDimensDummy( region,ilev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  mueT => region%levels(ilev)%turb%mueT
  tv   => region%levels(ilev)%mixt%tv
  gv   => region%levels(ilev)%mixt%gv
  tdv  => region%levels(ilev)%turb%dv

  DO k=kpcbeg-1,kpcend+1
    DO j=jpcbeg-1,jpcend+1
      DO i=ipcbeg-1,ipcend+1

        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        ijkNI = ijkN + 1
        ijkNJ = ijkN + iNOff
        ijkNK = ijkN + ijNOff

! ----- summing muet from the six faces has to be done here
!       since not all LES models calls TURB_LesCalcEddyVis
    
        tv(TV_MIXT_MUET,ijkC) = mueT(DIRI,ijkN)+mueT(DIRI,ijkNI) + &
                                mueT(DIRJ,ijkN)+mueT(DIRJ,ijkNJ) + &
                                mueT(DIRK,ijkN)+mueT(DIRK,ijkNK)
        tv(TV_MIXT_MUET,ijkC) = one6th*tv(TV_MIXT_MUET,ijkC)
        tdv(DV_TURB_CDYN,ijkC)= one6th*tdv(DV_TURB_CDYN,ijkC)

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif
#ifdef RFLU
  rPrt  =  1._RFREAL/region%mixtInput%prTurb
  indCp =  region%mixtInput%indCp
  mueT  => region%turb%mueT
  bMuet => region%turb%bMueT
  tv    => region%mixt%tv
  gv    => region%mixt%gv
  tdv   => region%turb%dv

! average muet at faces into cells 
! waiting for Andreas` routine for proper averaging

  tv(TV_MIXT_MUET,:) = 0._RFREAL

  DO ift = 1,region%grid%nFaces
    iC0 = region%grid%f2c(1,ift)
    iC1 = region%grid%f2c(2,ift)
    tv(TV_MIXT_MUET,iC0) = tv(TV_MIXT_MUET,iC0) + mueT(DIRI,ift)
    tv(TV_MIXT_MUET,iC1) = tv(TV_MIXT_MUET,iC1) + mueT(DIRI,ift)
  ENDDO    ! ift

  DO iPatch = 1,region%grid%nPatches
! TEMPORARY : removing usage of bf2bg from everywhere
!    ifgBeg = region%patches(iPatch)%bf2bg(BF2BG_BEG)
    ifg    = ifgBeg
    DO ifl = 1,region%patches(iPatch)%nBFaces
      iC0 = region%patches(iPatch)%bf2c(ifl)
      tv(TV_MIXT_MUET,iC0) = tv(TV_MIXT_MUET,iC0) + bMueT(DIRI,ifg)
      ifg = ifg + 1
    ENDDO    ! ifl
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

    tv( TV_MIXT_MUET,iC) = tv( TV_MIXT_MUET,iC)/nFacesPerCell
    tdv(DV_TURB_CDYN,iC) = tdv(DV_TURB_CDYN,iC)/nFacesPerCell
  ENDDO    ! iC
#endif

! extrapolate to dummy faces at the six sides of region; simple extrapolation 
! provides the original face TV values (through two point averaging) needed in 
! the computation of viscous fluxes

#ifdef RFLO
! extrapolate solution to dummy cells

  CALL TURB_FloExtrapolCellVec( region,TV_MIXT_MUET,TV_MIXT_MUET,tv )
#endif

! get turbulent thermal conductivity

#ifdef RFLO
  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend

        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
#endif
#ifdef RFLU
  DO iC = 1,region%grid%nCells
     ijkC = iC
#endif
        tv(TV_MIXT_MUET,ijkC) = MAX( tv(TV_MIXT_MUET,ijkC), 0._RFREAL )
        cpPrt = gv(GV_MIXT_CP,ijkC*indCp)*rPrt
        tv(TV_MIXT_TCOT,ijkC) = cpPrt*tv(TV_MIXT_MUET,ijkC)
#ifdef RFLO
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k
#endif
#ifdef RFLU
  ENDDO      ! iC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_GetTvCell

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_GetTvCell.F90,v $
! Revision 1.14  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2006/08/19 15:40:57  mparmar
! Removed bf2bg
!
! Revision 1.11  2006/01/13 03:47:11  wasistho
! initialize tv(muet,:) in Rocflu
!
! Revision 1.10  2005/12/30 23:21:29  wasistho
! bug fixed face indexing in patch
!
! Revision 1.9  2005/01/12 01:12:59  wasistho
! removed single quote signs since SUN has trouble with it
!
! Revision 1.8  2004/05/28 02:00:31  wasistho
! update unstructured grid LES
!
! Revision 1.7  2004/05/18 03:15:29  wasistho
! changed FloExtrapIntCellVec to FloExtrapolCellVec
!
! Revision 1.6  2004/05/17 20:47:48  wasistho
! compute first layer dummy mu_t instead of extrapolated
!
! Revision 1.5  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.4  2004/03/24 03:37:02  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/19 02:46:57  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.4  2004/02/18 02:23:11  wasistho
! used extrapolation routine to fill in dummy mut and klipped it
!
! Revision 1.3  2003/10/09 23:07:17  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!******************************************************************************








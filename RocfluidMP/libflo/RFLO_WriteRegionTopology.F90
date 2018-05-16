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
! Purpose: write topology of all regions to file.
!
! Description: none.
!
! Input: regions = region dimensions and topology
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_WriteRegionTopology.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_WriteRegionTopology( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(2*CHRLEN+4) :: fname

  INTEGER :: bcType, srcL1beg, srcL1end, srcL2beg, srcL2end, errorFlag

  TYPE(t_patch), POINTER :: patch

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_WriteRegionTopology',&
  'RFLO_WriteRegionTopology.F90' )

! open file & write number of regions

  fname = TRIM(global%outDir)//TRIM(global%casename)//'.top'
  OPEN(IF_TOPOL,file=fname,form='formatted',status='unknown',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  WRITE(IF_TOPOL,'(A)',err=10) '# topology for: '//TRIM(global%casename)
  WRITE(IF_TOPOL,'(A)',err=10) '#'
  WRITE(IF_TOPOL,  *  ,err=10) global%nRegions

! write topology of each region (at grid level 1 - the finest grid)

  DO iReg=1,global%nRegions
    WRITE(IF_TOPOL,*,err=10) iReg,regions(iReg)%nGridLevels
    WRITE(IF_TOPOL,*,err=10) regions(iReg)%nPatches, &
                                    regions(iReg)%levels(1)%grid%ipc, &
                                    regions(iReg)%levels(1)%grid%jpc, &
                                    regions(iReg)%levels(1)%grid%kpc

! - loop over all patches of a region

    DO iPatch=1,regions(iReg)%nPatches

      patch => regions(iReg)%levels(1)%patches(iPatch)
      bcType = patch%bcType

! --- patch with neighbour

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
          (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        IF (patch%align) THEN          ! l1 aligned with 11Src
          srcL1beg = -patch%srcL1beg
          srcL1end = -patch%srcL1end
          srcL2beg =  patch%srcL2beg
          srcL2end =  patch%srcL2end
        ELSE                           ! l1 aligned with l2Src
          srcL1beg =  patch%srcL1beg
          srcL1end =  patch%srcL1end
          srcL2beg = -patch%srcL2beg
          srcL2end = -patch%srcL2end
        ENDIF
        WRITE(IF_TOPOL,*,err=10)    &
          patch%bcType   ,patch%lbound   , &
          -patch%l1beg   ,-patch%l1end   , &
           patch%l2beg   , patch%l2end   , &
          patch%srcRegion,patch%srcLbound, &
          srcL1beg       ,srcL1end       , &
          srcL2beg       ,srcL2end       , &
          patch%bcCoupled

! --- no neighbour

      ELSE
        WRITE(IF_TOPOL,*,err=10)    &
          patch%bcType   ,patch%lbound   , &
          patch%l1beg    ,patch%l1end    , &
          patch%l2beg    ,patch%l2end    , &
          0,0,0,0,0,0    ,patch%bcCoupled
      ENDIF

    ENDDO   ! iPatch
  ENDDO     ! iReg

  CLOSE(IF_TOPOL,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! error handling --------------------------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

END SUBROUTINE RFLO_WriteRegionTopology

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_WriteRegionTopology.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.6  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.5  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.4  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.3  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
!******************************************************************************








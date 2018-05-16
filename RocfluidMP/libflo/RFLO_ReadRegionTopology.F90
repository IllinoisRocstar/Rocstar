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
! Purpose: read in topology of all regions (done on all processors);
!          figure out which l1,l2 directions are aligned on adjacent regions.
!
! Description: none.
!
! Input: none.
!
! Output: global%nRegions = number of regions
!         regions         = region dimensions and topology
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadRegionTopology.F90,v 1.12 2009/08/27 14:04:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadRegionTopology( global,regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(2*CHRLEN+4) :: fname
  CHARACTER(CHRLEN)     :: msg

  INTEGER :: regionNum, nLevels, nPatches, ipc, jpc, kpc, bcType, align
  INTEGER :: nCellsTot, errorFlag

  TYPE(t_patch), POINTER :: patch

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_ReadRegionTopology',&
  'RFLO_ReadRegionTopology.F90' )

! open file & read number of regions

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.top'
  OPEN(IF_TOPOL,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,&
    __LINE__,'File: '//TRIM(fname) )

  READ(IF_TOPOL,'(1X)',err=10,end=10)
  READ(IF_TOPOL,'(1X)',err=10,end=10)
  READ(IF_TOPOL,   *  ,err=10,end=10) global%nRegions

! allocate memory for region structure

  ALLOCATE( regions(global%nRegions),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! read topology of each region; store data at grid level 1 (finest grid)

  nCellsTot = 0

  DO iReg=1,global%nRegions
    READ(IF_TOPOL,*,err=10,end=10) regionNum,nLevels
    READ(IF_TOPOL,*,err=10,end=10) nPatches,ipc,jpc,kpc

    ALLOCATE( regions(regionNum)%levels(nLevels),stat=errorFlag )
    ALLOCATE( regions(regionNum)%levels(1)%patches(nPatches),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    regions(regionNum)%iRegionGlobal        = regionNum
    regions(regionNum)%nGridLevels          = nLevels
    regions(regionNum)%nPatches             = nPatches
    regions(regionNum)%levels(1)%grid%ipc   = ipc
    regions(regionNum)%levels(1)%grid%jpc   = jpc
    regions(regionNum)%levels(1)%grid%kpc   = kpc
    regions(regionNum)%active               = ACTIVE
    regions(regionNum)%mixtInput%externalBc = .false.

    nCellsTot = nCellsTot + ipc*jpc*kpc

! - set pointer "global" within regions

    regions(regionNum)%global => global

! - loop over all patches of a region

    DO iPatch=1,nPatches
      patch => regions(regionNum)%levels(1)%patches(iPatch)
      READ(IF_TOPOL,*,err=10,end=10)     &
        patch%bcType   ,patch%lbound   , &
        patch%l1beg    ,patch%l1end    , &
        patch%l2beg    ,patch%l2end    , &
        patch%srcRegion,patch%srcLbound, &
        patch%srcL1beg ,patch%srcL1end , &
        patch%srcL2beg ,patch%srcL2end , &
        patch%bcCoupled

      IF (patch%bcCoupled <= 0) THEN
        patch%bcCoupled = BC_INTERNAL
        patch%bcMotion  = BC_INTERNAL
      ELSE
        patch%bcCoupled = BC_EXTERNAL
        patch%bcMotion  = BC_EXTERNAL
      ENDIF
      IF (patch%bcType == BC_SLIPWALL_FREE   .OR. &
          patch%bcType == BC_SLIPWALL_FIXED  .OR. &
          patch%bcType == BC_SLIPWALL_XSLIDE .OR. &
          patch%bcType == BC_SLIPWALL_YSLIDE .OR. &
          patch%bcType == BC_SLIPWALL_ZSLIDE .OR. &
          patch%bcType == BC_SLIPWALL_XYSLIDE .OR. &
          patch%bcType == BC_SLIPWALL_XZSLIDE .OR. &
          patch%bcType == BC_SLIPWALL_YZSLIDE) THEN
        patch%bcMotion  = BC_EXTERNAL
      ENDIF
      IF (patch%bcType == BC_SYMMETRY_FREE   .OR. &
          patch%bcType == BC_SYMMETRY_FIXED  .OR. &
          patch%bcType == BC_SYMMETRY_XSLIDE .OR. &
          patch%bcType == BC_SYMMETRY_YSLIDE .OR. &
          patch%bcType == BC_SYMMETRY_ZSLIDE .OR. &
          patch%bcType == BC_SYMMETRY_XYSLIDE .OR. &
          patch%bcType == BC_SYMMETRY_XZSLIDE .OR. &
          patch%bcType == BC_SYMMETRY_YZSLIDE) THEN
        patch%bcMotion  = BC_EXTERNAL
      ENDIF
      IF (patch%bcType == BC_NOSLIPWALL_FREE   .OR. &
          patch%bcType == BC_NOSLIPWALL_FIXED  .OR. &
          patch%bcType == BC_NOSLIPWALL_XSLIDE .OR. &
          patch%bcType == BC_NOSLIPWALL_YSLIDE .OR. &
          patch%bcType == BC_NOSLIPWALL_ZSLIDE .OR. &
          patch%bcType == BC_NOSLIPWALL_XYSLIDE .OR. &
          patch%bcType == BC_NOSLIPWALL_XZSLIDE .OR. &
          patch%bcType == BC_NOSLIPWALL_YZSLIDE) THEN
        patch%bcMotion  = BC_EXTERNAL
      ENDIF
      IF (patch%bcType == BC_OUTFLOW_FREE   .OR. &
          patch%bcType == BC_OUTFLOW_FIXED  .OR. &
          patch%bcType == BC_OUTFLOW_XSLIDE .OR. &
          patch%bcType == BC_OUTFLOW_YSLIDE .OR. &
          patch%bcType == BC_OUTFLOW_ZSLIDE .OR. &
          patch%bcType == BC_OUTFLOW_XYSLIDE .OR. &
          patch%bcType == BC_OUTFLOW_XZSLIDE .OR. &
          patch%bcType == BC_OUTFLOW_YZSLIDE) THEN
        patch%bcMotion  = BC_EXTERNAL
      ENDIF
      IF (patch%bcCoupled == BC_EXTERNAL) &
        regions(regionNum)%mixtInput%externalBc = .true.

! --- check if BC type within range
      IF (patch%bcType<BC_CODE_MIN .OR. patch%bcType>BC_CODE_MAX) THEN
        WRITE(msg,'(A,I5,A)') 'Boundary code ',patch%bcType," ???"
        CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__,msg )
      ENDIF

! --- check if region face between 1 and 6
      IF (patch%lbound<1 .OR. patch%lbound>6) &
        CALL ErrorStop( global,ERR_WRONG_REGIONFACE,__LINE__ )

! --- initialize patch data
      patch%mixt%bcSet     = .true.          ! defaults: BC set
      patch%mixt%distrib   = BCDAT_CONSTANT  !   no distribution
      patch%mixt%nData     = 0               !   no data
      patch%mixt%nSwitches = 0               !   no switches

      patch%turb%bcSet     = .true.          ! defaults: BC set
      patch%turb%distrib   = BCDAT_CONSTANT  !   no distribution
      patch%turb%nData     = 0               !   no data
      patch%turb%nSwitches = 0               !   no switches

      patch%spec%bcSet     = .true.          ! defaults: BC set
      patch%spec%distrib   = BCDAT_CONSTANT  !   no distribution
      patch%spec%nData     = 0               !   no data
      patch%spec%nSwitches = 0               !   no switches

      patch%peul%bcSet     = .true.          ! defaults: BC set
      patch%peul%distrib   = BCDAT_CONSTANT  !   no distribution
      patch%peul%nData     = 0               !   no data
      patch%peul%nSwitches = 0               !   no switches

      patch%valRadi%bcSet     = .true.          ! defaults: BC set
      patch%valRadi%distrib   = BCDAT_CONSTANT  !   no distribution
      patch%valRadi%nData     = 0               !   no data
      patch%valRadi%nSwitches = 0               !   no switches

! --- unset BCs for certain types (where user input required)
      bcType = patch%bcType

      IF ((bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE) .OR. &
          (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
          (bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR. &
          (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE) .OR. &
          (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR. &
          (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE)) THEN
        patch%mixt%bcSet = .false.     ! will be set to true if BC defined
      ENDIF

      IF ((bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR. &
          (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR. &
          (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE)) THEN
        patch%peul%bcSet = .false.     ! will be set to true if BC defined
      ENDIF

    ENDDO   ! iPatch
  ENDDO     ! iReg

  CLOSE(IF_TOPOL,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! find aligned l1,l2 coordinates ----------------------------------------------

  DO iReg=1,global%nRegions
    DO iPatch=1,regions(iReg)%nPatches

      patch => regions(iReg)%levels(1)%patches(iPatch)

! --- patch with a neighbor

      bcType = patch%bcType
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
          (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN

! ----- check if source region within possible range

        IF (patch%srcRegion > global%nRegions) &
          CALL ErrorStop( global,ERR_REGION_RANGE,__LINE__ )

! ----- check number of source region face

        IF (patch%srcLbound<1 .OR. patch%srcLbound>6) &
          CALL ErrorStop( global,ERR_WRONG_REGIONFACE,__LINE__,'(source region).' )

! ----- present patch
        IF (patch%l1beg<0 .OR. patch%l1end<0 ) THEN        ! 1st direction
          align       = 10
          patch%l1beg = ABS(patch%l1beg)
          patch%l1end = ABS(patch%l1end)
          IF (patch%l2beg<0 .OR. patch%l2end<0 ) &
            CALL ErrorStop( global,ERR_PATCH_2ALIGN,__LINE__ )
        ELSE IF (patch%l2beg<0 .OR. patch%l2end<0 ) THEN   ! 2nd direction
          align       = 20
          patch%l2beg = ABS(patch%l2beg)
          patch%l2end = ABS(patch%l2end)
        ELSE
          CALL ErrorStop( global,ERR_PATCH_NOALIGN,__LINE__ )
        ENDIF

! ----- source patch
        IF (patch%srcL1beg<0 .OR. patch%srcL1end<0 ) THEN  ! 1st direction
          align          = align + 1
          patch%srcL1beg = ABS(patch%srcL1beg)
          patch%srcL1end = ABS(patch%srcL1end)
          IF (patch%srcL2beg<0 .OR. patch%srcL2end<0 ) &
            CALL ErrorStop( global,ERR_PATCH_2ALIGN,__LINE__ )
        ELSE IF (patch%srcL2beg<0 .OR. patch%srcL2end<0 ) THEN   ! 2nd dir.
          align          = align + 2
          patch%srcL2beg = ABS(patch%srcL2beg)
          patch%srcL2end = ABS(patch%srcL2end)
        ELSE
          CALL ErrorStop( global,ERR_PATCH_NOALIGN,__LINE__ )
        ENDIF

! ----- set alignment flag (1=yes, 0=no)
        IF (align==11 .OR. align==22) THEN
          patch%align = .true.
        ELSE
          patch%align = .false.
        ENDIF

! --- no neighbor

      ELSE
        patch%srcRegion = -999
        patch%srcPatch  = -999
        patch%align     = .true.
        patch%srcLbound = -999
        patch%srcL1beg  = -999
        patch%srcL1end  = -999
        patch%srcL2beg  = -999
        patch%srcL2end  = -999
      ENDIF

    ENDDO   ! iPatch
  ENDDO     ! iReg

! print some info -------------------------------------------------------------

  IF (global%myProcid==masterProc .AND. global%verbLevel>=VERBOSE_MED) THEN
    WRITE(STDOUT,'(/,A,I8)') SOLVER_NAME//' total no. of cells  = ',ncellsTot
    WRITE(STDOUT,'(A,I8,/)') SOLVER_NAME//' no. of grid regions = ', &
                             global%nRegions
  ENDIF

! error handling --------------------------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

END SUBROUTINE RFLO_ReadRegionTopology

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadRegionTopology.F90,v $
! Revision 1.12  2009/08/27 14:04:49  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.11  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/08/28 11:42:12  rfiedler
! Add grid motion constraint types for outflow BC.
!
! Revision 1.8  2006/08/24 13:15:36  rfiedler
! Rocflo now supports XYSLIDE, XZSLIDE, and YZSLIDE instead of TANGEN constraint.
!
! Revision 1.7  2006/08/19 15:38:22  mparmar
! Renamed patch variables
!
! Revision 1.6  2006/05/08 22:30:40  wasistho
! added prop-NS capability
!
! Revision 1.5  2006/05/04 04:24:28  wasistho
! include BC_SLIPWALL_FIXED in bcMotion externals
!
! Revision 1.4  2005/06/20 20:24:25  wasistho
! changed bcType to patch%bcType since bcType wasn't defined yet
!
! Revision 1.3  2005/06/19 05:31:16  wasistho
! shift index rocprop slipwalls
!
! Revision 1.2  2005/06/13 21:44:44  wasistho
! added new patch variable patch%bcMotion
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.14  2004/07/02 23:01:28  wasistho
! filled iRegionGlobal for Rocflo
!
! Revision 1.13  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.12  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.11  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.10  2002/09/27 03:34:04  jblazek
! BH seems to need "" around ??? in line 122.
!
! Revision 1.9  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.8  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.7  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.5  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.4  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.3  2002/03/18 21:56:39  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/27 18:38:19  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.1  2002/02/21 23:25:04  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.5  2002/01/11 17:13:30  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.4  2001/12/22 00:09:36  jblazek
! Added routines to store grid and solution.
!
! Revision 1.3  2001/12/19 23:09:20  jblazek
! Added routines to read grid and solution.
!
! Revision 1.2  2001/12/08 00:18:41  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************








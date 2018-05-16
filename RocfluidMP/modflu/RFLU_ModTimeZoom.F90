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
! ******************************************************************************
!
! Purpose: Collection of routines related to solution zooming..
!
! Description: None.
!
! Notes: None. 
!**********************************************************************

MODULE RFLU_ModTimeZoom

  USE ModDataTypes   !add use list
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI

  IMPLICIT NONE

  TYPE t_TimeZoom
     REAL(RFREAL),Dimension(CV_MIXT_DENS:CV_MIXT_ENER) :: CvBulk,resbar_Qbulk,&
          CvdVdtBulk
  END TYPE t_TimeZoom
  
  PRIVATE
  PUBLIC :: RFLU_TimeZoomDriver,t_TimeZoom,TimeZoom, &
       RFLU_UnZoomGridSpeeds,RFLU_ZoomGridSpeeds


!Module Local varaibles
  TYPE (t_TimeZoom) :: TimeZoom

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  CHARACTER(CHRLEN) :: & 
       RCSIdentString = '$RCSfile: RFLU_ModTimeZoom.F90,v $' 

! ******************************************************************************
! Routines
! ******************************************************************************

CONTAINS


  SUBROUTINE RFLU_TimeZoomDriver(regions)
    IMPLICIT NONE

!---Dummy
    TYPE(t_region), POINTER :: regions(:)

    CALL RFLU_TimeZoomComputeBulkVars(regions)
    CALL RFLU_TimeZoomSumAddSource(regions)

    RETURN
  END SUBROUTINE RFLU_TimeZoomDriver


  SUBROUTINE RFLU_TimeZoomComputeBulkVars(regions)
    IMPLICIT NONE

!---Dummy
    TYPE(t_region), POINTER :: regions(:)
!---Local
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_global), POINTER :: global
    REAL(RFREAL), DIMENSION(:), POINTER :: vol,volOld
    REAL(RFREAL), DIMENSION(:,:), POINTER :: Cv
    INTEGER :: iReg,iCell,i,nVars,nCvvars
    REAL(RFREAL) ::x, LocalIntegVol,GlobalIntegVol,dvdt_i, volBulkOld
    REAL(RFREAL),DIMENSION(2*(CV_MIXT_ENER-CV_MIXT_DENS+1)+2) :: BulkVarsLocal,BulkVarsGlobal
!---For submerged nozzles
    REAL(RFREAL) :: y, z
    REAL(RFREAL) :: RNozzleInlet2

    global => regions(1)%global

    RNozzleInlet2 = global%tzNozRad**2

    nVars = UBOUND(BulkVarsLocal,1)-LBOUND(BulkVarsLocal,1)+1
    nCvvars = CV_MIXT_ENER-CV_MIXT_DENS+1
    BulkVarsLocal = 0_RFREAL

    DO iReg = 1,global%nRegionsLocal

! Region Pointers
       pRegion => regions(iReg)
       vol    => pRegion%grid%vol
       volOld => pRegion%gridOld%vol
       Cv     => pRegion%mixt%cv

! Sum over Cells
       IF (.NOT. global%tzSubNoz) THEN
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)
             IF(x > global%tzMinPlane .AND. x < global%tzMaxPlane) THEN
                dvdt_i = (vol(iCell)-volOld(iCell))/MAX(global%dtMin,  &
                         TINY(1.0_RFREAL))
                BulkVarsLocal  = BulkVarsLocal + (/ (Cv(i,iCell)*vol(iCell), &
                     i=CV_MIXT_DENS,CV_MIXT_ENER),vol(iCell),volOld(iCell), &
                     (Cv(i,iCell)*dvdt_i,i=CV_MIXT_DENS,CV_MIXT_ENER)/)
             ENDIF
          END DO !iCell
       ELSE
!         Exclude submerged nozzle
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)
             y = pRegion%grid%cofg(global%tzCoordTrans1,iCell)
             z = pRegion%grid%cofg(global%tzCoordTrans2,iCell)
             IF((x > global%tzMinPlane .AND. x < global%tzNozInlet) .OR.  &
                ((x >= global%tzNozInlet .AND. x < global%tzMaxPlane) .AND.  &
                 (y**2 + z**2 > RNozzleInlet2) ) ) THEN
                dvdt_i = (vol(iCell)-volOld(iCell))/MAX(global%dtMin,  &
                         TINY(1.0_RFREAL))
                BulkVarsLocal  = BulkVarsLocal + (/ (Cv(i,iCell)*vol(iCell), &
                     i=CV_MIXT_DENS,CV_MIXT_ENER),vol(iCell),volOld(iCell), &
                     (Cv(i,iCell)*dvdt_i,i=CV_MIXT_DENS,CV_MIXT_ENER)/)
             ENDIF ! In zoomed volume
          END DO !iCell
       ENDIF ! submerged
    END DO  !iReg

    CALL MPI_ALLREDUCE(BulkVarsLocal, BulkVarsGlobal,nVars, &
	         MPI_RFREAL,MPI_SUM,global%mpiComm,global%error)
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    
    global%tzVolBulk = BulkVarsGlobal(nCvvars+1)
    volBulkOld = BulkVarsGlobal(nCvvars+2)
    TimeZoom%CvBulk(CV_MIXT_DENS:CV_MIXT_ENER) = BulkVarsGlobal(1:nCvvars)/global%tzVolBulk
    TimeZoom%CvdVdtBulk(CV_MIXT_DENS:CV_MIXT_ENER) = BulkVarsGlobal(nCvvars+3:nVars)

    global%tzRadChamber = sqrt(global%tzVolBulk/global%tzLenChamb/global%pi)
    global%tzEpsNozz     = (global%tzThroatRad/global%tzRadChamber)**2.0D0
   
!   is the time difference between pRegion%gridOld and pRegion%grid = global%dtMin !???ZOOM???
    global%tzDvolBulkDt = (global%tzVolBulk-volBulkOld)/max(global%dtMin,TINY(1.0_RFREAL))


    return
  End Subroutine RFLU_TimeZoomComputeBulkVars



  Subroutine RFLU_TimeZoomSumAddSource(regions)
    USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed,&
         &                            RFLU_ScaleGridSpeed

    IMPLICIT NONE

!---Dummy
    TYPE(t_region), POINTER :: regions(:)
!---Local
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: pRegion
    REAL(RFREAL), DIMENSION(:,:), POINTER :: res   !use res to conform to Q1D code
    REAL(RFREAL), DIMENSION(:), POINTER :: vol,volOld
    INTEGER :: iReg,iCell,nVars,ivar
    REAL(RFREAL) :: x,betaFill,volBulk,ddt_volBulk,termR1,termR2,termR3,termL1
!---For submerged nozzles
    REAL(RFREAL) :: y, z
    REAL(RFREAL) :: RNozzleInlet2
    
    global => regions(1)%global

    RNozzleInlet2 = global%tzNozRad**2

    nVars  = CV_MIXT_ENER-CV_MIXT_DENS+1

    betaFill    = global%zoomFactor
    volBulk     = global%tzVolBulk
    ddt_volBulk = global%tzDvolBulkDt


!sum Residuals
    call RFLU_TimeZoomSumResiduals(regions)

!Add R1 
    DO iReg = 1,global%nRegionsLocal
       pRegion => regions(iReg)
       res     => pRegion%mixt%rhs
       vol     => pRegion%grid%vol
       volOld  => pRegion%gridOld%vol

       IF (.NOT. global%tzSubNoz) THEN
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)

             if(x > global%tzMinPlane .and. x < global%tzMaxPlane) then
                termR1 = (1_RFREAL-betaFill)/betaFill/max(global%dtMin,TINY(1.0_RFREAL))
                termR1 = RFLU_ScaleGridSpeed(pRegion,termR1)
                do ivar = CV_MIXT_DENS,CV_MIXT_ENER
                   res(ivar,iCell) = res(ivar,iCell) + termR1*pRegion%mixt%Cv(ivar,iCell)*&
                        &(vol(iCell)-volOld(iCell))
                end do
             end if ! if x > Xmin & x < Xmax

          End DO !ic
       ELSE
!         Exclude submerged nozzle
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)
             y = pRegion%grid%cofg(global%tzCoordTrans1,iCell)
             z = pRegion%grid%cofg(global%tzCoordTrans2,iCell)
             IF((x > global%tzMinPlane .AND. x < global%tzNozInlet) .OR.  &
                ((x >= global%tzNozInlet .AND. x < global%tzMaxPlane) .AND.  &
                 (y**2 + z**2 > RNozzleInlet2) ) ) THEN

                termR1 = (1_RFREAL-betaFill)/betaFill/max(global%dtMin,TINY(1.0_RFREAL))
                termR1 = RFLU_ScaleGridSpeed(pRegion,termR1)
                do ivar = CV_MIXT_DENS,CV_MIXT_ENER
                   res(ivar,iCell) = res(ivar,iCell) + termR1*pRegion%mixt%Cv(ivar,iCell)*&
                        &(vol(iCell)-volOld(iCell))
                end do
             ENDIF ! In zoomed volume

          End DO !ic
       ENDIF ! submerged
    End Do !iReg

! Add R2 and R3
    DO iReg = 1,global%nRegionsLocal
       pRegion => regions(iReg)
       res    => pRegion%mixt%rhs
       vol    => pRegion%grid%vol

       IF (.NOT. global%tzSubNoz) THEN
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)

             if(x > global%tzMinPlane .and. x < global%tzMaxPlane) then
   !---R2
                termR2 = (betaFill-1_RFREAL)*vol(iCell)/volBulk
                do ivar = CV_MIXT_DENS,CV_MIXT_ENER
                   res(ivar,iCell) = res(ivar,iCell) + termR2*TimeZoom%resbar_Qbulk(ivar)
                end do
   !---R3
                termR3 = (betaFill-1_RFREAL)/betaFill*vol(iCell)/volBulk*ddt_volBulk 
                termR3 = RFLU_ScaleGridSpeed(pRegion,termR3)
                do ivar = CV_MIXT_DENS,CV_MIXT_ENER
                   res(ivar,iCell) = res(ivar,iCell) + termR3*TimeZoom%CvBulk(ivar)
                end do

             end if ! if x > Xmin & x < Xmax

          End DO !ic
       ELSE
!         Exclude submerged nozzle
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)
             y = pRegion%grid%cofg(global%tzCoordTrans1,iCell)
             z = pRegion%grid%cofg(global%tzCoordTrans2,iCell)
             IF((x > global%tzMinPlane .AND. x < global%tzNozInlet) .OR.  &
                ((x >= global%tzNozInlet .AND. x < global%tzMaxPlane) .AND.  &
                 (y**2 + z**2 > RNozzleInlet2) ) ) THEN

!---R2
                termR2 = (betaFill-1_RFREAL)*vol(iCell)/volBulk
                do ivar = CV_MIXT_DENS,CV_MIXT_ENER
                   res(ivar,iCell) = res(ivar,iCell) + termR2*TimeZoom%resbar_Qbulk(ivar)
                end do
!---R3
                termR3 = (betaFill-1_RFREAL)/betaFill*vol(iCell)/volBulk*ddt_volBulk 
                termR3 = RFLU_ScaleGridSpeed(pRegion,termR3)
                do ivar = CV_MIXT_DENS,CV_MIXT_ENER
                   res(ivar,iCell) = res(ivar,iCell) + termR3*TimeZoom%CvBulk(ivar)
                end do
             ENDIF ! In zoomed volume

          End DO !ic

       ENDIF ! submerged

    End Do !iReg

    return
  End Subroutine RFLU_TimeZoomSumAddSource
    



  Subroutine RFLU_TimeZoomSumResiduals(Regions)
    IMPLICIT NONE

!---Dummy
    TYPE(t_global), POINTER :: global
!---Local
    TYPE(t_region), POINTER :: pRegion, regions(:) 
    REAL(RFREAL), DIMENSION(:,:), POINTER :: res   !use res to conform to Q1D code
    INTEGER :: iReg,iCell,nVars
    REAL(RFREAL) ::x, LocalIntegVol,GlobalIntegVol
    REAL(RFREAL),DIMENSION(CV_MIXT_ENER-CV_MIXT_DENS+1) :: ResBarVecLocal,ResBarVecGlobal
!---For submerged nozzles
    REAL(RFREAL) :: y, z
    REAL(RFREAL) :: RNozzleInlet2

    global => regions(1)%global

    RNozzleInlet2 = global%tzNozRad**2

    nVars = CV_MIXT_ENER-CV_MIXT_DENS+1

    ResBarVecLocal = 0_RFREAL
    DO iReg = 1,global%nRegionsLocal

       pRegion => regions(iReg)
       res    => pRegion%mixt%rhs

       IF (.NOT. global%tzSubNoz) THEN
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)
             if(x > global%tzMinPlane .and. x < global%tzMaxPlane) then
                ResBarVecLocal = ResBarVecLocal + res(CV_MIXT_DENS:CV_MIXT_ENER,iCell)
             end if
          End DO !iCell
       ELSE 
!         Exclude submerged nozzle
          DO iCell = 1,pRegion%grid%nCells
             x = pRegion%grid%cofg(global%tzCoordLong,iCell)
             y = pRegion%grid%cofg(global%tzCoordTrans1,iCell)
             z = pRegion%grid%cofg(global%tzCoordTrans2,iCell)
             IF((x > global%tzMinPlane .AND. x < global%tzNozInlet) .OR.  &
                ((x >= global%tzNozInlet .AND. x < global%tzMaxPlane) .AND.  &
                 (y**2 + z**2 > RNozzleInlet2) ) ) THEN
                ResBarVecLocal = ResBarVecLocal + res(CV_MIXT_DENS:CV_MIXT_ENER,iCell)
             ENDIF ! In zoomed volume
          End DO !iCell
       ENDIF ! submerged
    End DO

    CALL MPI_ALLREDUCE(ResBarVecLocal, ResBarVecGlobal,nVars, &
         MPI_RFREAL,MPI_SUM,global%mpiComm,global%error)
    IF ( global%error /= ERR_NONE ) THEN
       CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    
    TimeZoom%resbar_Qbulk(CV_MIXT_DENS:CV_MIXT_ENER) = ResBarVecGlobal(1:nvars)



    return
  END SUBROUTINE RFLU_TimeZoomSumResiduals



! ******************************************************************************
!
! Purpose: Time UnZoom grid speeds for faces and boundary patches.
!
! Description: None.
!
! Input:
!   pRegion      	Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_UnZoomGridSpeeds(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: ifg,ifl,iPatch
    REAL(RFREAL) :: scaleFactor
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! *****************************************************************************
!   Set pointers and variables
! *****************************************************************************

    global => pRegion%global
    pGrid => pRegion%grid

    scaleFactor = 1_RFREAL/global%Zoomfactor
  
! ******************************************************************************
!   Scale grid speeds
! ******************************************************************************

    DO ifg = LBOUND(pGrid%gs,1),UBOUND(pGrid%gs,1) 
      pGrid%gs(ifg) = scaleFactor*pGrid%gs(ifg)
    END DO ! ifg
 
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifl = LBOUND(pPatch%gs,1),UBOUND(pPatch%gs,1) 
        pPatch%gs(ifl) = scaleFactor*pPatch%gs(ifl)
      END DO ! ifl
    END DO ! iPatch
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_UnZoomGridSpeeds




! ******************************************************************************
!
! Purpose: Time Zoom grid speeds for faces and boundary patches.
!
! Description: None.
!
! Input:
!   pRegion      	Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ZoomGridSpeeds(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: ifg,ifl,iPatch
    REAL(RFREAL) :: scaleFactor
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! *****************************************************************************
!   Set pointers and variables
! *****************************************************************************

    global => pRegion%global
    pGrid => pRegion%grid

    scaleFactor = global%Zoomfactor
  
! ******************************************************************************
!   Scale grid speeds
! ******************************************************************************

    DO ifg = LBOUND(pGrid%gs,1),UBOUND(pGrid%gs,1) 
      pGrid%gs(ifg) = scaleFactor*pGrid%gs(ifg)
    END DO ! ifg
 
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifl = LBOUND(pPatch%gs,1),UBOUND(pPatch%gs,1) 
        pPatch%gs(ifl) = scaleFactor*pPatch%gs(ifl)
      END DO ! ifl
    END DO ! iPatch
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ZoomGridSpeeds



! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModTimeZoom


! ******************************************************************************
!
! RCS Revision history:
!
!
! ******************************************************************************
  







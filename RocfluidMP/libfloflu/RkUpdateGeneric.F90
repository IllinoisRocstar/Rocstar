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
! Purpose: Updates solution with classical 4-stage Runge-Kutta method.
!
! Description: None.
!
! Input: 
!   region      Region data
!   varType     Variable type to be updated
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   cv          Conserved variables
!   cvOld       Old conserved variables
!   rhs         Residual
!   rhsSum      Residual sum
!
! Output: 
!   cv          Conserved variables
!   rhsSum      Residual sum
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkUpdateGeneric.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkUpdateGeneric(region,varType,iStage,icBeg,icEnd,ivBeg,ivEnd, &
                           cv,cvOld,rhs,rhsSum)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region
  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd,varType
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,rhs,rhsSum

! =============================================================================
! Locals
! =============================================================================

  LOGICAL :: moveGrid
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
#ifdef RFLO
  INTEGER :: iLev
#endif  
  REAL(RFREAL) :: adtv,fac,volRat  
  REAL(RFREAL) :: ark(5),grk(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: vol,volOld
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RkUpdateGeneric.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction(global,'RkUpdateGeneric',&
  'RkUpdateGeneric.F90')

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! =============================================================================
! Set volume(s) when updating cell-based variables
! =============================================================================

  IF ( varType == VAR_TYPE_CELL ) THEN 
    moveGrid = region%mixtInput%moveGrid

#ifdef RFLO
    iLev = region%currLevel

    vol => region%levels(iLev)%grid%vol

    IF ( moveGrid .EQV. .TRUE. ) THEN 
      volOld => region%levels(iLev)%gridOld%vol
    END IF ! moveGrid
#endif

#ifdef RFLU
    vol => region%grid%vol

    IF ( moveGrid .EQV. .TRUE. ) THEN 
      volOld => region%gridOld%vol
    END IF ! moveGrid
#endif
  END IF ! varType

! *****************************************************************************
! Update
! *****************************************************************************

  fac = ark(iStage)*global%dtMin

  SELECT CASE ( varType ) 
  
! =============================================================================
!   Update cell-based variables, for which we need the volume (and volume 
!   ratio for moving grid computations). 
! =============================================================================
  
    CASE ( VAR_TYPE_CELL )   
        
! -----------------------------------------------------------------------------
!     Update for moving grid
! -----------------------------------------------------------------------------

      IF ( moveGrid .EQV. .TRUE. ) THEN
        SELECT CASE ( global%rkScheme ) 
          CASE ( RK_SCHEME_4_CLASSICAL ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic      
            ELSE IF ( iStage == global%nrkSteps ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd                                       
                  cv(iv,ic) = volRat*cvOld(iv,ic) & 
                            - adtv*(rhs(iv,ic)+rhsSum(iv,ic))
                END DO ! iv
              END DO ! ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhsSum(iv,ic) + grk(iStage)*rhs(iv,ic)
                END DO ! iv
              END DO ! ic
            END IF ! iStage
          CASE ( RK_SCHEME_3_WRAY ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic      
            ELSE IF ( iStage == 2 ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd                                       
                  cv(iv,ic)     = volRat*cvOld(iv,ic) & 
                                - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                  rhsSum(iv,ic) = rhs(iv,ic)              
                END DO ! iv
              END DO ! ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic) = volRat*cvOld(iv,ic) & 
                            - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                END DO ! iv
              END DO ! ic
            END IF ! iStage      
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! global%rkScheme

! =============================================================================
!     Update for non-moving grid
! =============================================================================

      ELSE 
        SELECT CASE ( global%rkScheme ) 
          CASE ( RK_SCHEME_4_CLASSICAL ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd                               
                  cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic
            ELSE IF ( iStage == global%nrkSteps ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic) = cvOld(iv,ic) - adtv*(rhs(iv,ic) + rhsSum(iv,ic))
                END DO ! iv
              END DO !ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd          
                  cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhsSum(iv,ic) + grk(iStage)*rhs(iv,ic)
                END DO ! iv
              END DO !ic      
            END IF ! iStage
          CASE ( RK_SCHEME_3_WRAY ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd                               
                  cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic
            ELSE IF ( iStage == 2 ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = cvOld(iv,ic) & 
                                - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                  rhsSum(iv,ic) = rhs(iv,ic)              
                END DO ! iv
              END DO !ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd          
                  cv(iv,ic) = cvOld(iv,ic) & 
                            - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                END DO ! iv
              END DO !ic      
            END IF ! iStage      
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! global%rkScheme
      END IF ! moveGrid

! =============================================================================
!   Update point-based variables, for which we DO NOT need the volume (and 
!   volume ratio for moving grid computations). 
! =============================================================================
  
    CASE ( VAR_TYPE_POINT )
      adtv = fac

      SELECT CASE ( global%rkScheme ) 
        CASE ( RK_SCHEME_4_CLASSICAL ) 
          IF ( iStage == 1 ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd                               
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                rhsSum(iv,ic) = rhs(iv,ic)
              END DO ! iv
            END DO ! ic
          ELSE IF ( iStage == global%nrkSteps ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd
                cv(iv,ic) = cvOld(iv,ic) - adtv*(rhs(iv,ic) + rhsSum(iv,ic))
              END DO ! iv
            END DO !ic      
          ELSE
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd          
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                rhsSum(iv,ic) = rhsSum(iv,ic) + grk(iStage)*rhs(iv,ic)
              END DO ! iv
            END DO !ic      
          END IF ! iStage
        CASE ( RK_SCHEME_3_WRAY ) 
          IF ( iStage == 1 ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd                               
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                rhsSum(iv,ic) = rhs(iv,ic)
              END DO ! iv
            END DO ! ic
          ELSE IF ( iStage == 2 ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd
                cv(iv,ic)     = cvOld(iv,ic) & 
                              - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                rhsSum(iv,ic) = rhs(iv,ic)              
              END DO ! iv
            END DO !ic      
          ELSE
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd          
                cv(iv,ic) = cvOld(iv,ic) & 
                          - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
              END DO ! iv
            END DO !ic      
          END IF ! iStage      
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! global%rkScheme
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
  END SELECT ! varType

! *****************************************************************************
! End
! *****************************************************************************

 CALL DeregisterFunction(global)

END SUBROUTINE RkUpdateGeneric

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkUpdateGeneric.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:51:08  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/11/17 16:24:30  haselbac
! Added varType and RK3
!
! Revision 1.1  2003/11/25 21:01:50  haselbac
! Initial revision
!
!******************************************************************************








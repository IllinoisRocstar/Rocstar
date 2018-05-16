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
! Purpose: Compute integrals 1,2,4,5 of optimal LES approach.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to current region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_AllocateDCUHREArrays.F90,v 1.5 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************


    SUBROUTINE RFLU_AllocateDCUHREArrays(nDim,nFunNZ,nFun)
    
      IMPLICIT NONE
      
! --- parameters      
      
      INTEGER, INTENT(IN) :: nDim,nFunNZ,nFun        

! --- locals

      INTEGER :: errorFlag

! ==============================================================================
!     Start
! ==============================================================================
      
      CALL RegisterFunction( 'RFLU_AllocateDCUHREArrays',&
  'RFLU_AllocateDCUHREArrays.F90' )

! ==============================================================================
!     Allocate memory
! ==============================================================================

      ALLOCATE(lowLim(nDim),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_ALLOCATE,__LINE__,'lowLim')
      END IF ! global%error

      ALLOCATE(uppLim(nDim),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_ALLOCATE,__LINE__,'uppLim')
      END IF ! global%error

      ALLOCATE(errAbsEst(nFun),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_ALLOCATE,__LINE__,'errAbsEst')
      END IF ! global%error

      ALLOCATE(integralNZ(nFunNZ),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_ALLOCATE,__LINE__,'integralNZ')
      END IF ! global%error

! ==============================================================================
!     End
! ==============================================================================

      CALL DeregisterFunction        
    
    END SUBROUTINE RFLU_AllocateDCUHREArrays





! ******************************************************************************
!   Deallocate DCUHRE arrays
! ******************************************************************************

    SUBROUTINE RFLU_DeallocateDCUHREArrays
    
      IMPLICIT NONE
      
! --- locals      
      
      INTEGER :: errorFlag
      
! ==============================================================================
!     Start
! ==============================================================================
      
      CALL RegisterFunction( 'RFLU_DeallocateDCUHREArrays',&
  'RFLU_AllocateDCUHREArrays.F90' )

! ==============================================================================
!     Allocate memory
! ==============================================================================

      DEALLOCATE(lowLim,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_DEALLOCATE,__LINE__,'lowLim')
      END IF ! global%error

      DEALLOCATE(uppLim,STAT=errorFlag)
      global%error = errorFlag      
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_DEALLOCATE,__LINE__,'uppLim')
      END IF ! global%error

      DEALLOCATE(errAbsEst,STAT=errorFlag)
      global%error = errorFlag      
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_DEALLOCATE,__LINE__,'errAbsEst')
      END IF ! global%error

      DEALLOCATE(integralNZ,STAT=errorFlag)
      global%error = errorFlag      
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(ERR_DEALLOCATE,__LINE__,'integralNZ')
      END IF ! global%error

! ==============================================================================
!     End
! ==============================================================================

      CALL DeregisterFunction        
    
    END SUBROUTINE RFLU_DeallocateDCUHREArrays









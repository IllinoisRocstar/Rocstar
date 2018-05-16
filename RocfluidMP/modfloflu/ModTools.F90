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
! Purpose: Collection of utility functions.
!
! Description: None
!
! Notes: 
!   1. Named ModTools instead of ModUtilities because we already have a 
!      utilities directory.
!
!******************************************************************************
!
! $Id: ModTools.F90,v 1.11 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
!******************************************************************************

MODULE ModTools

  USE ModDataTypes

  IMPLICIT NONE
  
  CONTAINS
  
! ------------------------------------------------------------------------------
!   Swap two integers
! ------------------------------------------------------------------------------

    SUBROUTINE SwapIntegers(a,b)
   
      INTEGER, INTENT(INOUT) :: a,b
 
      INTEGER :: c
    
      c = b
      b = a
      a = c
      
    END SUBROUTINE SwapIntegers
    
! ------------------------------------------------------------------------------
!   Swap two floats 
! ------------------------------------------------------------------------------

    SUBROUTINE SwapRFREALs(a,b)
   
      REAL(RFREAL), INTENT(INOUT) :: a,b
      
      REAL(RFREAL) :: c
    
      c = b
      b = a
      a = c
    
    END SUBROUTINE SwapRFREALs
  
! ------------------------------------------------------------------------------
!   Prevent division by zero by blending of variable with machine precision
! ------------------------------------------------------------------------------

    FUNCTION MakeNonZero(x)
   
      REAL(RFREAL) :: MakeNonZero
   
      REAL(RFREAL), INTENT(IN) :: x
    
      MakeNonZero = x + SIGN(1.0_RFREAL,x)*EPSILON(1.0_RFREAL)
    
    END FUNCTION MakeNonZero

! ------------------------------------------------------------------------------
!   Comparison of floating-point numbers for equality
! ------------------------------------------------------------------------------

    LOGICAL FUNCTION FloatEqual(a,b,tolIn)
   
      REAL(RFREAL), INTENT(IN) :: a,b
      REAL(RFREAL), INTENT(IN), OPTIONAL :: tolIn
   
      REAL(RFREAL) :: tol 
    
      floatEqual = .FALSE.      

      IF ( PRESENT(tolIn) .EQV. .TRUE. ) THEN 
        tol = tolIn
      ELSE  
        tol = 10.0_RFREAL*EPSILON(1.0_RFREAL)
      END IF ! PRESENT

      IF ( ABS(a-b) <= (1.0_RFREAL + 0.5_RFREAL*ABS(a+b))*tol ) THEN 
        floatEqual = .TRUE.
      END IF ! ABS
        
    END FUNCTION FloatEqual

! ------------------------------------------------------------------------------
!   Comparison of floating-point numbers: greater than
! ------------------------------------------------------------------------------

    LOGICAL FUNCTION FloatGreater(a,b)
   
      REAL(RFREAL), INTENT(IN) :: a,b
    
      floatGreater = .FALSE.
    
      IF ( a - b > b*EPSILON(1.0_RFREAL) ) THEN 
        floatGreater = .TRUE.
      END IF ! a
        
    END FUNCTION FloatGreater

! ------------------------------------------------------------------------------
!   Comparison of floating-point numbers: less than
! ------------------------------------------------------------------------------

    LOGICAL FUNCTION FloatLess(a,b)
   
      REAL(RFREAL), INTENT(IN) :: a,b
    
      floatLess = .FALSE.
    
      IF ( a - b < b*EPSILON(1.0_RFREAL) ) THEN 
        floatLess = .TRUE.
      END IF ! a
        
    END FUNCTION FloatLess

! ------------------------------------------------------------------------------
!   Compute Factorial
! ------------------------------------------------------------------------------

    INTEGER FUNCTION CompFact(n)
    
      INTEGER, INTENT(IN) :: n
      
      INTEGER :: i
    
      CompFact = 1

      DO i = 2,n
        CompFact = CompFact*i
      END DO ! i    
    
    END FUNCTION CompFact 

! ------------------------------------------------------------------------------
!   Compute (a^2 + b^2)**(1/2) without under- or overflow
! ------------------------------------------------------------------------------

    FUNCTION CompPythag(a,b)
     
      REAL(RFREAL) :: CompPythag
    
      REAL(RFREAL), INTENT(IN) :: a,b
      
      REAL(RFREAL) :: absa,absb
      
      absa = ABS(a)
      absb = ABS(b)

      IF ( absa > absb ) THEN 
        CompPythag = absa*SQRT(1.0_RFREAL + (absb/absa)**2)
      ELSE
        IF ( absb == 0.0_RFREAL ) THEN 
          CompPythag = 0.0_RFREAL
        ELSE 
          CompPythag = absb*SQRT(1.0_RFREAL + (absa/absb)**2)
        END IF ! absb  
      END IF ! absa      
    
    END FUNCTION CompPythag

! ------------------------------------------------------------------------------
! Detect NaNs
! ------------------------------------------------------------------------------

  FUNCTION IsNan(x)
  
    LOGICAL :: IsNan
    
    REAL(RFREAL), INTENT(IN) :: x
    
    IsNan = .FALSE. 
    
    IF ( .NOT.(x > -HUGE(1.0) .AND. x < HUGE(1.0)) ) THEN
      IsNan = .TRUE.
    END IF ! NOT
  
  END FUNCTION IsNan

! ******************************************************************************
! End
! ******************************************************************************

END MODULE ModTools

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModTools.F90,v $
! Revision 1.11  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/03/25 21:47:54  haselbac
! Renamed SwapFloats to SwapRFREALS
!
! Revision 1.8  2003/12/04 03:28:29  haselbac
! Added function to detect NaNs
!
! Revision 1.7  2003/02/06 19:30:59  haselbac
! Added optional tolerance argument to FloatEqual
!
! Revision 1.6  2003/01/28 16:49:01  haselbac
! Changed FloatEqual to work better with very small floats
!
! Revision 1.5  2002/11/27 20:24:42  haselbac
! Changed tolerance in FloatEqual
!
! Revision 1.4  2002/11/26 15:25:45  haselbac
! Fixed bug in FloatEqual
!
! Revision 1.3  2002/09/09 15:00:54  haselbac
! Added robust calculation of CompPythag
!
! Revision 1.2  2002/07/25 15:13:39  haselbac
! Added factorial function
!
! Revision 1.1  2002/05/04 16:51:58  haselbac
! Initial revision
!
!******************************************************************************







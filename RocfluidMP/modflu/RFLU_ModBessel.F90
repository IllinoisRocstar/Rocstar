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
! Purpose: Routines to compute Bessel functions, their derivatives, and their
!   zeros.
!
! Description: None.
!
! Notes: 
!   1. From: http://perso.wanadoo.fr/jean-pierre.moreau/f_bessel.html, 
!      labelled: F90 Release 1.0 By J-P Moreau, Paris. The original program
!      header included the following reference:
!      Ref.: www.esrf.fr/computing/expg/libraries/smf/PROGRAMS/MJYZO.FOR
!
! ******************************************************************************
!
! $Id: RFLU_ModBessel.F90,v 1.4 2008/12/06 08:44:20 mtcampbe Exp $
!
! ******************************************************************************

MODULE RFLU_ModBessel
      
  USE ModDataTypes    
      
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_JYNDD, & 
            RFLU_JYZO, & 
            RFLU_JYZOM

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModBessel.F90,v $ $Revision: 1.4 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose:  Compute Bessel functions Jn(x) and Yn(x) and their first and 
!   second derivatives 
!
! Description: None.
!
! Input:
!   n           Order of Jn(x) and Yn(x)
!   x           Argument of Jn(x) and Yn(x) (x>0)
!
! Output:
!   BJN         Jn(x)
!   DJN         dJn/dx(x)
!   FJN         d2Jn/dx2(x)
!   BYN         Yn(x)
!   DYN         dYn/dx(x)
!   FYN         d2Yn/dx2(x)
!
! Notes:
!   1. Removed labelled loops and CONTINUE statements from original code. 
!
! ******************************************************************************
        
  SUBROUTINE RFLU_JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: N
    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION, INTENT(OUT) :: BJN,DJN,FJN,BYN,DYN,FYN

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: K,M,MT,NT
    DOUBLE PRECISION :: BS,EC,E0,F,F0,F1,SU,S1
    DOUBLE PRECISION :: BJ(102),BY(102)

! ******************************************************************************
!   Start
! ******************************************************************************
                
    DO NT=1,900
       MT=INT(0.5*LOG10(6.28*NT)-NT*LOG10(1.36*DABS(X)/NT))
       IF (MT.GT.20) THEN 
         EXIT
       END IF ! MT
    END DO ! NT    
    M=NT
    BS=0.0D0
    F0=0.0D0
    F1=1.0D-35
    SU=0.0D0    
    DO K=M,0,-1
      F=2.0D0*(K+1.0D0)*F1/X-F0       
      IF (K.LE.N+1) THEN 
        BJ(K+1)=F
      END IF ! K       
      IF (K.EQ.2*INT(K/2)) THEN
         BS=BS+2.0D0*F
         IF (K.NE.0) THEN
           SU=SU+(-1)**(K/2)*F/K
         END IF ! K
      END IF       
      F0=F1 
      F1=F
    END DO ! K    
    DO K=0,N+1
      BJ(K+1)=BJ(K+1)/(BS-F)
    END DO ! K 
    BJN=BJ(N+1)
    EC=0.5772156649015329D0
    E0=0.3183098861837907D0
    S1=2.0D0*E0*(DLOG(X/2.0D0)+EC)*BJ(1)
    F0=S1-8.0D0*E0*SU/(BS-F)
    F1=(BJ(2)*F0-2.0D0*E0/X)/BJ(1)
    BY(1)=F0
    BY(2)=F1    
    DO K=2,N+1
       F=2.0D0*(K-1.0D0)*F1/X-F0
       BY(K+1)=F
       F0=F1
       F1=F
    END DO ! K    
    BYN=BY(N+1)
    DJN=-BJ(N+2)+N*BJ(N+1)/X
    DYN=-BY(N+2)+N*BY(N+1)/X
    FJN=(N*N/(X*X)-1.0D0)*BJN-DJN/X
    FYN=(N*N/(X*X)-1.0D0)*BYN-DYN/X

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_JYNDD




! ******************************************************************************
!
! Purpose: Compute the zeros of Bessel functions Jn(x) and Yn(x) and their 
!   first derivatives
!
! Description: None.
!
! Input:
!   N           Order of Bessel functions (0 to 100)
!   NT          Number of zeros
!
! Output:
!   RJ0(L)      L-th zero of Jn(x), L=1,2,...,NT
!   RJ1(L)      L-th zero of dJn/dx(x), L=1,2,...,NT
!   RY0(L)      L-th zero of Yn(x), L=1,2,...,NT
!   RY1(L)      L-th zero of dYn/dx(x), L=1,2,...,NT
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JYZO(N,NT,RJ0,RJ1,RY0,RY1)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: N,NT
    DOUBLE PRECISION, INTENT(OUT) :: RJ0(NT),RJ1(NT),RY0(NT),RY1(NT)    

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: L
    DOUBLE PRECISION :: BJN,BYN,DJN,DYN,FJN,FYN,X,X0

! ******************************************************************************
!   Start
! ******************************************************************************
        
    IF (N.LE.20) THEN
       X=2.82141+1.15859*N
    ELSE
       X=N+1.85576*N**0.33333+1.03315/N**0.33333
    ENDIF
    L=0
10  X0=X
    CALL RFLU_JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
    X=X-BJN/DJN
    IF (DABS(X-X0).GT.1.0D-9) GO TO 10
    L=L+1
    RJ0(L)=X
    X=X+3.1416+(0.0972+0.0679*N-0.000354*N**2)/L
    IF (L.LT.NT) GO TO 10
    IF (N.LE.20) THEN
       X=0.961587+1.07703*N
    ELSE
       X=N+0.80861*N**0.33333+0.07249/N**0.33333
    ENDIF
    IF (N.EQ.0) X=3.8317
    L=0
15  X0=X
    CALL RFLU_JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
    X=X-DJN/FJN
    IF (DABS(X-X0).GT.1.0D-9) GO TO 15
    L=L+1
    RJ1(L)=X
    X=X+3.1416+(0.4955+0.0915*N-0.000435*N**2)/L
    IF (L.LT.NT) GO TO 15
    IF (N.LE.20) THEN
       X=1.19477+1.08933*N
    ELSE
       X=N+0.93158*N**0.33333+0.26035/N**0.33333
    ENDIF           
    L=0
20  X0=X
    CALL RFLU_JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
    X=X-BYN/DYN
    IF (DABS(X-X0).GT.1.0D-9) GO TO 20
    L=L+1
    RY0(L)=X
    X=X+3.1416+(0.312+0.0852*N-0.000403*N**2)/L
    IF (L.LT.NT) GO TO 20
    IF (N.LE.20) THEN
       X=2.67257+1.16099*N
    ELSE
       X=N+1.8211*N**0.33333+0.94001/N**0.33333
    ENDIF  
    L=0
25  X0=X
    CALL RFLU_JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
    X=X-DYN/FYN
    IF (DABS(X-X0).GT.1.0D-9) GO TO 25
    L=L+1
    RY1(L)=X
    X=X+3.1416+(0.197+0.0643*N-0.000286*N**2)/L 
    IF (L.LT.NT) GO TO 25

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_JYZO
        
        
        
        
        
! ******************************************************************************
!
! Purpose: Compute Mth zero of Bessel functions Jn(x) and Yn(x) and their 
!   first derivatives
!
! Description: None.
!
! Input:
!   N           Order of Bessel functions (0 to 100)
!   M           Index of zero
!
! Output:
!   RJ0M        M-th zero of Jn(x)
!   RJ1M        M-th zero of dJn/dx(x)
!   RY0M        M-th zero of Yn(x)
!   RY1M        M-th zero of dYn/dx(x)
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JYZOM(N,M,RJ0M,RJ1M,RY0M,RY1M)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: N,M
    DOUBLE PRECISION, INTENT(OUT) :: RJ0M,RJ1M,RY0M,RY1M    

! ==============================================================================
!   Locals
! ==============================================================================

    DOUBLE PRECISION :: BJN,BYN,DJN,DYN,FJN,FYN,X,X0
    DOUBLE PRECISION :: RJ0(M),RJ1(M),RY0(M),RY1(M) 

! ******************************************************************************
!   Start
! ******************************************************************************
        
    CALL RFLU_JYZO(N,M,RJ0,RJ1,RY0,RY1)
        
    RJ0M = RJ0(M)
    RJ1M = RJ1(M)
    RY0M = RY0(M)
    RY1M = RY1(M)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_JYZOM       
        
     


! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModBessel


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModBessel.F90,v $
! Revision 1.4  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.1  2005/03/15 20:42:53  haselbac
! Initial revision
!
! ******************************************************************************







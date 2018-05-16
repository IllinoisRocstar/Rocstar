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
!*******************************************************************************
!
! Purpose: Suite of routines to carry out octree operations.
!
! Description: None.
!
! Notes: 
!   1. The routines contained in this module are from Tim Baker. They are
!      modified so that they no longer use COMMON statements and follow the 
!      design philosophy of Rocflu. Some of the routines have been renamed:
!
!      The routine RFLU_BuildOctree was called OCT_SETUP by Tim Baker, and 
!      RFLU_QueryOctree was called OCT_SEARCH.   
!
!      To create and use a octree, one has to take the following steps:
!        1. RFLU_CreateOctree
!        2. RFLU_BuildOctree
!        3. RFLU_QueryOctree
!        4. RFLU_DestroyOctree
!
!*******************************************************************************
!
! $Id: RFLU_ModOctree.F90,v 1.9 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModOctree

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_CreateOctree,RFLU_BuildOctree,RFLU_QueryOctree, & 
            RFLU_DestroyOctree
  
  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
  INTEGER, PARAMETER :: MTEST =  5000, & 
                        MNUMP =   100
  
  INTEGER :: MOCTR,NNODE
       
  INTEGER, DIMENSION(:), ALLOCATABLE :: ICHK,NLINK
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: NOCTR      
       
  REAL(RFREAL) :: DISMIN     
  REAL(RFREAL), DIMENSION(2) :: XFAR,YFAR,ZFAR 
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: X,Y,Z    
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
! ==============================================================================
!   Create octree
! ==============================================================================  

    SUBROUTINE RFLU_CreateOctree(global,nPoints)
    
      IMPLICIT NONE   
        
      INTEGER, INTENT(IN) :: nPoints 
      TYPE(t_global), POINTER :: global
        
      INTEGER :: errorFlag  
        
      CALL RegisterFunction(global,'RFLU_CreateOctree',&
  'RFLU_ModOctree.F90')

      NNODE = nPoints
      MOCTR = 2*NNODE

      ALLOCATE(ICHK(NNODE),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ICHK')
      END IF ! global%error     

      ALLOCATE(NLINK(NNODE),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'NLINK')
      END IF ! global%error

      ALLOCATE(NOCTR(2,MOCTR),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'NOCTR')
      END IF ! global%error
      
      ALLOCATE(X(NNODE),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'X')
      END IF ! global%error

      ALLOCATE(Y(NNODE),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'Y')
      END IF ! global%error

      ALLOCATE(Z(NNODE),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'Z')
      END IF ! global%error
    
      CALL DeregisterFunction(global)    

    END SUBROUTINE RFLU_CreateOctree


! ==============================================================================
!   Build octree
! ==============================================================================  

      SUBROUTINE RFLU_BuildOctree(XI,YI,ZI,XLOW,XUPP,YLOW,YUPP,ZLOW,ZUPP)
!
!     ******************************************************************
!     *                                                                *
!     *  GIVEN A MESH DEFINED AT THE SET OF POINTS                     *
!     *  (X(N),Y(N),Z(N), N=1,NNODE), STORE THEM IN AN OCTREE DATA     *
!     *  STRUCTURE FOR USE IN FUTURE SEARCH QUERIES.                   *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   2002                            *
!     *                                                                *
!     ******************************************************************
!
!     ******************************************************************
!

! --- arguments

      REAL(RFREAL), INTENT(IN) :: XLOW,XUPP,YLOW,YUPP,ZLOW,ZUPP
      REAL(RFREAL), INTENT(IN) :: XI(NNODE),YI(NNODE),ZI(NNODE)

! --- locals

      INTEGER :: N
!
!     ******************************************************************
!
! --- copy coordinates into internal coordinate array
!
      DO N=1,NNODE      
        X(N) = XI(N)
        Y(N) = YI(N)
        Z(N) = ZI(N)
      END DO ! N      
!
!     DETERMINE MINIMUM AND MAXIMUM EXTENT OF ORIGINAL DOMAIN
!
      ICHK(1:NNODE) = 0
      
      XFAR(1) = XLOW
      XFAR(2) = XUPP
      YFAR(1) = YLOW
      YFAR(2) = YUPP
      ZFAR(1) = ZLOW
      ZFAR(2) = ZUPP            
!
!     INSERT ORIGINAL MESH POINTS INTO OCTREE DATA STRUCTURE
!
      NLINK(1)   = 0
      NOCTR(1,1) = 1
      NOCTR(2,1) = 0
      
      CALL OCTREE   
      
      END SUBROUTINE RFLU_BuildOctree



! =============================================================================
!     Construct octree
! =============================================================================

      SUBROUTINE OCTREE
!
!      ******************************************************************
!      *                                                                *
!      *   SET UP OCTREE DATA STRUCTURE FOR THE POINTS (X(N),Y(N),Z(N)) *
!      *   WHERE N RANGES FROM 1 TO NNODE                               *
!      *                                                                *
!      ******************************************************************
!      ******************************************************************
!      *                                                                *
!      *   COPYRIGHT (C) TIM BAKER     2002                             *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!
! --- locals 

      INTEGER :: I,IOCTR,IROOT,J,K,L,LTEST,N,NEXT,NTEST,NXSGN,NYSGN,NZSGN
      INTEGER :: JJ(8),NSTORE(8)
      REAL(RFREAL) :: TOL, & 
                      XHALF,XHIGH,XLOW,XSHIFT,XSIZE,XSGN, & 
                      YHALF,YHIGH,YLOW,YSHIFT,YSIZE,YSGN, & 
                      ZHALF,ZHIGH,ZLOW,ZSHIFT,ZSIZE,ZSGN
      REAL(RFREAL) :: XROOT(2),YROOT(2),ZROOT(2)
!
!     ******************************************************************
! 
      TOL        = 1.000001_RFREAL
      IOCTR      = 1
!
!     LOOP OVER LIST OF POINTS
!
      DO 150 N=2,NNODE
      I          = 1
      XLOW       = XFAR(1)
      XHIGH      = XFAR(2)
      YLOW       = YFAR(1)
      YHIGH      = YFAR(2)
      ZLOW       = ZFAR(1)
      ZHIGH      = ZFAR(2)    
      
      IF (NOCTR(1,I).LT.0) GO TO 80
      L          = 0
      NEXT       = NOCTR(1,I) 
    6 NEXT       = NLINK(NEXT)
      L          = L  +1  
      IF (NEXT.NE.0) GO TO 6
   10 IROOT      = I
      XROOT(1)   = XLOW
      XROOT(2)   = XHIGH
      YROOT(1)   = YLOW
      YROOT(2)   = YHIGH
      ZROOT(1)   = ZLOW
      ZROOT(2)   = ZHIGH
      IF (L.EQ.8) GO TO 30
      L          = L  +1
      NLINK(N)   = 0
      IF (L.GT.1) GO TO 15
      NOCTR(1,I) = N
      GO TO 150
   15 NEXT       = NOCTR(1,I)
   20 NTEST      = NLINK(NEXT)
      IF (NTEST.EQ.0) GO TO 25
      NEXT       = NTEST
      GO TO 20
   25 NLINK(NEXT) = N
      GO TO 150
!
!     OCTANT NOW CONTAINS NINE POINTS. FORM EIGHT NEW OCTANTS
!     AND AMEND TREE STRUCTURE.
!
   30 DO 35 K=1,8
      IOCTR      = IOCTR  +1
      IF (IOCTR.GT.MOCTR) GO TO 210
      NOCTR(2,IOCTR) = I
   35 NOCTR(1,IOCTR) = 0
!
!     ASSIGN POINTS TO NEW OCTANTS
!
      NEXT       = NOCTR(1,I)
      DO 40 L=1,8
      NSTORE(L)  = NEXT
      NEXT       = NLINK(NEXT)      
   40 JJ(L)      = 0
      DO 45 K=1,8
      J          = NSTORE(K)
      XSHIFT     = X(J)  -.5_RFREAL*(XLOW  +XHIGH)
      XSIZE      = MAX(1.E-9_RFREAL,ABS(XSHIFT))
      NXSGN      = (INT(TOL*XSHIFT/XSIZE)  +1)/2*2
      YSHIFT     = Y(J)  -.5_RFREAL*(YLOW  +YHIGH)
      YSIZE      = MAX(1.E-9_RFREAL,ABS(YSHIFT))
      NYSGN      = (INT(TOL*YSHIFT/YSIZE)  +1)/2*2
      ZSHIFT     = Z(J)  -.5_RFREAL*(ZLOW  +ZHIGH)
      ZSIZE      = MAX(1.E-9_RFREAL,ABS(ZSHIFT))
      NZSGN      = (INT(TOL*ZSHIFT/ZSIZE)  +1)/2*2
      L          = 1  +NXSGN/2  +NYSGN  +2*NZSGN
      JJ(L)      = JJ(L)  +1
      NLINK(J)   = 0
      IF (JJ(L).GT.1) GO TO 41
      NOCTR(1,IOCTR-8+L) = J
      GO TO 45
   41 NEXT       = NOCTR(1,IOCTR-8+L)
   42 NTEST      = NLINK(NEXT)
      IF (NTEST.EQ.0) GO TO 43
      NEXT       = NTEST
      GO TO 42
   43 NLINK(NEXT) = J
   45 CONTINUE
!
!     CHECK WHETHER ALL EIGHT POINTS LIE IN ONE OCTANT
!
      NOCTR(1,I) = 7  -IOCTR
      LTEST      = 0
      DO 60 L=1,8
      IF (JJ(L).EQ.8) LTEST = L
   60 CONTINUE
      IF (LTEST.EQ.0) GO TO 70
      I          = IOCTR  -8  +LTEST
      XHALF      = .5_RFREAL*(XLOW  +XHIGH)
      XSGN       = MOD(LTEST,2)
      XLOW       = XSGN*XLOW   +(1._RFREAL  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1._RFREAL  -XSGN)*XHIGH
      YHALF      = .5_RFREAL*(YLOW  +YHIGH)
      YSGN       = ISIGN(1,2*MOD(LTEST-1,4)-3)
      YLOW       = .5_RFREAL*((1._RFREAL-YSGN)*YLOW +(1._RFREAL+YSGN)*YHALF)
      YHIGH      = .5_RFREAL*((1._RFREAL-YSGN)*YHALF+(1._RFREAL+YSGN)*YHIGH)
      ZHALF      = .5_RFREAL*(ZLOW  +ZHIGH)
      ZSGN       = ISIGN(1,2*LTEST-9)
      ZLOW       = .5_RFREAL*((1._RFREAL-ZSGN)*ZLOW +(1._RFREAL+ZSGN)*ZHALF)
      ZHIGH      = .5_RFREAL*((1._RFREAL-ZSGN)*ZHALF+(1._RFREAL+ZSGN)*ZHIGH)
      GO TO 30
!
!     LOCATE SUB-OCTANT IN WHICH POINT LIES
!
   70 I          = IROOT
      XLOW       = XROOT(1)
      XHIGH      = XROOT(2)
      YLOW       = YROOT(1)
      YHIGH      = YROOT(2)
      ZLOW       = ZROOT(1)
      ZHIGH      = ZROOT(2)
   80 XHALF      = .5_RFREAL*(XLOW  +XHIGH)
      YHALF      = .5_RFREAL*(YLOW  +YHIGH)
      ZHALF      = .5_RFREAL*(ZLOW  +ZHIGH)
      XSHIFT     = X(N)  -XHALF
      XSIZE      = MAX(1.E-9_RFREAL,ABS(XSHIFT))
      NXSGN      = (INT(TOL*XSHIFT/XSIZE)  +1)/2*2
      YSHIFT     = Y(N)  -YHALF
      YSIZE      = MAX(1.E-9_RFREAL,ABS(YSHIFT))
      NYSGN      = (INT(TOL*YSHIFT/YSIZE)  +1)/2*2
      ZSHIFT     = Z(N)  -ZHALF
      ZSIZE      = MAX(1.E-9_RFREAL,ABS(ZSHIFT))
      NZSGN      = (INT(TOL*ZSHIFT/ZSIZE)  +1)/2*2
      L          = 1  +NXSGN/2  +NYSGN  +2*NZSGN
      I          = -NOCTR(1,I)  +L  -1
      XSGN       = MOD(L,2)
      XLOW       = XSGN*XLOW   +(1._RFREAL  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1._RFREAL  -XSGN)*XHIGH
      YSGN       = ISIGN(1,2*MOD(L-1,4)-3)
      YLOW       = .5_RFREAL*((1._RFREAL-YSGN)*YLOW +(1._RFREAL+YSGN)*YHALF)
      YHIGH      = .5_RFREAL*((1._RFREAL-YSGN)*YHALF+(1._RFREAL+YSGN)*YHIGH)
      ZSGN       = ISIGN(1,2*L-9)
      ZLOW       = .5_RFREAL*((1._RFREAL-ZSGN)*ZLOW +(1._RFREAL+ZSGN)*ZHALF)
      ZHIGH      = .5_RFREAL*((1._RFREAL-ZSGN)*ZHALF+(1._RFREAL+ZSGN)*ZHIGH)
      IF (NOCTR(1,I).LT.0) GO TO 80
      L          = 0
      IF (NOCTR(1,I).EQ.0) GO TO 10
      NEXT       = NOCTR(1,I)
   85 NEXT       = NLINK(NEXT)
      L          = L  +1
      IF (NEXT.NE.0) GO TO 85
      GO TO 10
  150 CONTINUE
      RETURN
  210 WRITE (6,610)
  
      STOP

  610 FORMAT(5X,'DIMENSION OF NOCTR ARRAY EXCEEDED IN ROUTINE OCTREE.'/ &
             5X,'INCREASE SIZE OF MOCTR')
     
      END SUBROUTINE OCTREE


! ==============================================================================
!   Query octree
! ==============================================================================

      SUBROUTINE RFLU_QueryOctree(XPT,YPT,ZPT,NUMP,NEIGHP)
!
!     ******************************************************************
!     *                                                                *
!     *  FIND THE NUMP POINTS IN THE OCTREE THAT ARE CLOSEST TO THE    *
!     *  POINT (XPT,YPT,ZPT). THE ADDRESSES OF THE CLOSEST POINTS ARE  *
!     *  GIVEN BY (NEIGHP(M), M = 1, NUMP).                            *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *   COPYRIGHT (C) TIM BAKER     2002                             *
!     *                                                                *
!     ******************************************************************
!
!     ******************************************************************
!

! --- arguments

      INTEGER, INTENT(IN) :: NUMP
      INTEGER, INTENT(INOUT) :: NEIGHP(NUMP)
      REAL(RFREAL), INTENT(IN) :: XPT,YPT,ZPT

! --- locals

      INTEGER :: I,IC,IFLAG,II,ITRY,ISTART,J,JJ,K,KC,L,LFLAG,M, & 
                 NXSGN,NYSGN,NZSGN
      INTEGER :: NSRCH(MTEST),KSRCH(MTEST)
      REAL(RFREAL) :: DIST,DMIN,TOL,TOLPT, &
                      XHALF,XHIGH,XHIGH0,XL,XLOW,XLOW0,XSHIFT,XSGN,XSIZE,XU, & 
                      YHALF,YHIGH,YHIGH0,YL,YLOW,YLOW0,YSHIFT,YSGN,YSIZE,YU, & 
                      ZHALF,ZHIGH,ZHIGH0,ZL,ZLOW,ZLOW0,ZSHIFT,ZSGN,ZSIZE,ZU
      REAL(RFREAL) ::  XOCTR(2,MTEST),YOCTR(2,MTEST),ZOCTR(2,MTEST), &
                       XHOLD(2,MTEST),YHOLD(2,MTEST),ZHOLD(2,MTEST), &
                       XKEEP(2),YKEEP(2),ZKEEP(2)
!
!     ******************************************************************
!
      TOLPT      = 1.E-9_RFREAL
      TOL        = 1.000001_RFREAL      
      IF (XPT.LT.XFAR(1)-TOLPT.OR.XPT.GT.XFAR(2)+TOLPT) GO TO 200
      IF (YPT.LT.YFAR(1)-TOLPT.OR.YPT.GT.YFAR(2)+TOLPT) GO TO 200
      IF (ZPT.LT.ZFAR(1)-TOLPT.OR.ZPT.GT.ZFAR(2)+TOLPT) GO TO 200
!
!     FIRST FIND OCTANT WHICH CONTAINS POINT (XPT,YPT,ZPT)
!
      M          = 0
      I          = 1
      XLOW       = XFAR(1)
      XHIGH      = XFAR(2)
      YLOW       = YFAR(1)
      YHIGH      = YFAR(2)
      ZLOW       = ZFAR(1)
      ZHIGH      = ZFAR(2)
      XLOW0      = XLOW
      XHIGH0     = XHIGH
      YLOW0      = YLOW
      YHIGH0     = YHIGH
      ZLOW0      = ZLOW
      ZHIGH0     = ZHIGH      
      ISTART     = 1
      IF (NOCTR(1,I).GT.0) GO TO 50
!
!     LOCATE SUB-OCTANT IN WHICH POINT LIES
!
   20 XHALF      = .5_RFREAL*(XLOW  +XHIGH)
      YHALF      = .5_RFREAL*(YLOW  +YHIGH)
      ZHALF      = .5_RFREAL*(ZLOW  +ZHIGH)
      XSHIFT     = XPT  -XHALF
      XSIZE      = MAX(1.E-9_RFREAL,ABS(XSHIFT))
      NXSGN      = (INT(TOL*XSHIFT/XSIZE)  +1)/2*2
      YSHIFT     = YPT  -YHALF
      YSIZE      = MAX(1.E-9_RFREAL,ABS(YSHIFT))
      NYSGN      = (INT(TOL*YSHIFT/YSIZE)  +1)/2*2
      ZSHIFT     = ZPT  -ZHALF
      ZSIZE      = MAX(1.E-9_RFREAL,ABS(ZSHIFT))
      NZSGN      = (INT(TOL*ZSHIFT/ZSIZE)  +1)/2*2
      L          = 1  +NXSGN/2  +NYSGN  +2*NZSGN
      I          = -NOCTR(1,I)  +L  -1
      XSGN       = MOD(L,2)
      XLOW       = XSGN*XLOW   +(1.  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1.  -XSGN)*XHIGH
      YSGN       = ISIGN(1,2*MOD(L-1,4)-3)
      YLOW       = .5_RFREAL*((1._RFREAL-YSGN)*YLOW +(1._RFREAL+YSGN)*YHALF)
      YHIGH      = .5_RFREAL*((1._RFREAL-YSGN)*YHALF+(1._RFREAL+YSGN)*YHIGH)
      ZSGN       = ISIGN(1,2*L-9)
      ZLOW       = .5_RFREAL*((1._RFREAL-ZSGN)*ZLOW +(1._RFREAL+ZSGN)*ZHALF)
      ZHIGH      = .5_RFREAL*((1._RFREAL-ZSGN)*ZHALF+(1._RFREAL+ZSGN)*ZHIGH)
      IF (NOCTR(1,I).LT.0) GO TO 20
      XLOW0      = XLOW
      XHIGH0     = XHIGH
      YLOW0      = YLOW
      YHIGH0     = YHIGH
      ZLOW0      = ZLOW
      ZHIGH0     = ZHIGH
      ISTART     = I
!
!     SEARCH FOR POINT NEIGHP(M) THAT IS NEAREST TO (XPT,YPT,ZPT)
!
   50 M           = M  +1  
      NEIGHP(M)   = 1
      DISMIN      = (XFAR(2)  -XFAR(1))**2  +(YFAR(2)  -YFAR(1))**2 &
                   +(ZFAR(2)  -ZFAR(1))**2
      XLOW       = XLOW0
      XHIGH      = XHIGH0
      YLOW       = YLOW0
      YHIGH      = YHIGH0
      ZLOW       = ZLOW0
      ZHIGH      = ZHIGH0
      I          = ISTART
      IF (NOCTR(1,I).EQ.0) GO TO 65
      J           = NOCTR(1,I)
   60 DIST        = (XPT  -X(J))**2  +(YPT  -Y(J))**2 &
                   +(ZPT  -Z(J))**2
      IF (DIST.LT.DISMIN.AND.ICHK(J).EQ.0) THEN
         NEIGHP(M)   = J
         DISMIN      = DIST
      ENDIF
      J           = NLINK(J)
      IF (J.NE.0) GO TO 60
   65 DMIN        = SQRT(DISMIN)
      XL          = XPT  -DMIN
      XU          = XPT  +DMIN
      YL          = YPT  -DMIN
      YU          = YPT  +DMIN
      ZL          = ZPT  -DMIN
      ZU          = ZPT  +DMIN
   70 IF (XL.GT.XLOW.AND.XU.LT.XHIGH.AND. &
          YL.GT.YLOW.AND.YU.LT.YHIGH.AND. &
          ZL.GT.ZLOW.AND.ZU.LT.ZHIGH) THEN
         IF (M.EQ.NUMP) THEN
            DO 72 M=1,NUMP
            ICHK(NEIGHP(M)) = 0
   72       CONTINUE
            RETURN
         ELSE
            ICHK(NEIGHP(M)) = 1
            GO TO 50
         ENDIF
      ENDIF
      IF (I.EQ.1) THEN
         IF (M.EQ.NUMP) THEN
            DO 74 M=1,NUMP
            ICHK(NEIGHP(M)) = 0
   74       CONTINUE
            RETURN
         ELSE
            ICHK(NEIGHP(M)) = 1
            GO TO 50
         ENDIF
      ENDIF
      IFLAG       = I
      I           = NOCTR(2,I)
      LFLAG       = IFLAG  +NOCTR(1,I)  +1
      XSGN        = MOD(LFLAG,2)
      XOCTR(1,1)  = (2._RFREAL  -XSGN)*XLOW  -(1._RFREAL  -XSGN)*XHIGH
      XOCTR(2,1)  = (1._RFREAL  +XSGN)*XHIGH  -XSGN*XLOW
      YSGN        = ISIGN(1,2*MOD(LFLAG-1,4)-3)
      YOCTR(1,1)  = .5_RFREAL*((3._RFREAL+YSGN)*YLOW  -(1._RFREAL+YSGN)*YHIGH)
      YOCTR(2,1)  = .5_RFREAL*((3._RFREAL-YSGN)*YHIGH -(1._RFREAL-YSGN)*YLOW)
      ZSGN        = ISIGN(1,2*LFLAG-9)
      ZOCTR(1,1)  = .5_RFREAL*((3._RFREAL+ZSGN)*ZLOW  -(1._RFREAL+ZSGN)*ZHIGH)
      ZOCTR(2,1)  = .5_RFREAL*((3._RFREAL-ZSGN)*ZHIGH -(1._RFREAL-ZSGN)*ZLOW)
      XKEEP(1)    = XOCTR(1,1)
      XKEEP(2)    = XOCTR(2,1)
      YKEEP(1)    = YOCTR(1,1)
      YKEEP(2)    = YOCTR(2,1)
      ZKEEP(1)    = ZOCTR(1,1)
      ZKEEP(2)    = ZOCTR(2,1)
!
!     EXAMINE PARENT OF OCTANT AND ALL ITS OFFSPRING
!
      KC          = 1
      KSRCH(1)    = I
   75 IC          = 0
      DO 90 J=1,KC
      II          = KSRCH(J)
      DO 90 K=1,8
      ITRY        = -NOCTR(1,II)  +K  -1
      IF (ITRY.EQ.IFLAG) GO TO 90
      XHALF      = .5_RFREAL*(XOCTR(1,J)  +XOCTR(2,J))
      XSGN       = MOD(K,2)
      XLOW       = XSGN*XOCTR(1,J)   +(1._RFREAL  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1._RFREAL  -XSGN)*XOCTR(2,J)
      YHALF      = .5_RFREAL*(YOCTR(1,J)  +YOCTR(2,J))
      YSGN       = ISIGN(1,2*MOD(K-1,4)-3)
      YLOW       = .5_RFREAL*((1._RFREAL-YSGN)*YOCTR(1,J) +(1._RFREAL+YSGN)*YHALF)
      YHIGH      = .5_RFREAL*((1._RFREAL-YSGN)*YHALF  +(1._RFREAL+YSGN)*YOCTR(2,J))
      ZHALF      = .5_RFREAL*(ZOCTR(1,J)  +ZOCTR(2,J))
      ZSGN       = ISIGN(1,2*K-9)
      ZLOW       = .5_RFREAL*((1._RFREAL-ZSGN)*ZOCTR(1,J)+(1._RFREAL+ZSGN)*ZHALF)
      ZHIGH      = .5_RFREAL*((1._RFREAL-ZSGN)*ZHALF +(1._RFREAL+ZSGN)*ZOCTR(2,J))
      IF (XL.GT.XHIGH.OR.XU.LT.XLOW.OR. &
          YL.GT.YHIGH.OR.YU.LT.YLOW.OR. &
          ZL.GT.ZHIGH.OR.ZU.LT.ZLOW) GO TO 90
      IF (NOCTR(1,ITRY).GE.0) GO TO 80
      IC          = IC  +1
      IF (IC.GT.MTEST) GO TO 210
      XHOLD(1,IC) = XLOW
      XHOLD(2,IC) = XHIGH
      YHOLD(1,IC) = YLOW
      YHOLD(2,IC) = YHIGH
      ZHOLD(1,IC) = ZLOW
      ZHOLD(2,IC) = ZHIGH
      NSRCH(IC)   = ITRY
      GO TO 90
   80 IF (NOCTR(1,ITRY).EQ.0) GO TO 90
      JJ          = NOCTR(1,ITRY)
   85 DIST        = (XPT  -X(JJ))**2  +(YPT  -Y(JJ))**2 &
                   +(ZPT  -Z(JJ))**2
      IF (DIST.LT.DISMIN.AND.ICHK(JJ).EQ.0) THEN
         NEIGHP(M)   = JJ
         DISMIN      = DIST
      ENDIF
      JJ          = NLINK(JJ)
      IF (JJ.NE.0) GO TO 85
      DMIN        = SQRT(DISMIN)
      XL          = XPT  -DMIN
      XU          = XPT  +DMIN
      YL          = YPT  -DMIN
      YU          = YPT  +DMIN
      ZL          = ZPT  -DMIN
      ZU          = ZPT  +DMIN
   90 CONTINUE
      IF (IC.EQ.0) THEN
         XLOW        = XKEEP(1)
         XHIGH       = XKEEP(2)
         YLOW        = YKEEP(1)
         YHIGH       = YKEEP(2)
         ZLOW        = ZKEEP(1)
         ZHIGH       = ZKEEP(2)
         GO TO 70
      ENDIF
      KC          = IC
      DO 95 J=1,KC
      KSRCH(J)    = NSRCH(J)
      XOCTR(1,J)  = XHOLD(1,J)
      XOCTR(2,J)  = XHOLD(2,J)
      YOCTR(1,J)  = YHOLD(1,J)
      YOCTR(2,J)  = YHOLD(2,J)
      ZOCTR(1,J)  = ZHOLD(1,J)
      ZOCTR(2,J)  = ZHOLD(2,J)
   95 CONTINUE
      GO TO 75
  200 WRITE (6,600) XPT,YPT,ZPT
      STOP
  210 WRITE (6,610)
      STOP
      
  600 FORMAT(5X,'POINT XPT = ',F12.4,' YPT = ',F12.4,' ZPT = ',F12.4/ &
             5X,'LIES OUTSIDE CONVEX HULL'/ &
             5X,'PROGRAM STOPPED IN ROUTINE OCT_SEARCH') 
  610 FORMAT(5X,'DIMENSION OF NSRCH AND KSRCH ARRAYS EXCEEDED'/ &
             5X,'IN ROUTINE OCT_SEARCH. INCREASE SIZE OF MTEST')
     
      END SUBROUTINE RFLU_QueryOctree


! ==============================================================================
!   Destroy octree
! ==============================================================================  
   
    SUBROUTINE RFLU_DestroyOctree(global)
    
      IMPLICIT NONE

      TYPE(t_global), POINTER :: global

      INTEGER :: errorFlag

      CALL RegisterFunction(global,'RFLU_DestroyOctree',&
  'RFLU_ModOctree.F90')

      DEALLOCATE(ICHK,STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ICHK')
      END IF ! global%error     

      DEALLOCATE(NLINK,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'NLINK')
      END IF ! global%error

      DEALLOCATE(NOCTR,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'NOCTR')
      END IF ! global%error
 
      DEALLOCATE(X,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'X')
      END IF ! global%error

      DEALLOCATE(Y,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'Y')
      END IF ! global%error

      DEALLOCATE(Z,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'Z')
      END IF ! global%error
 
      CALL DeregisterFunction(global)

    END SUBROUTINE RFLU_DestroyOctree


 
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModOctree


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModOctree.F90,v $
!   Revision 1.9  2008/12/06 08:44:23  mtcampbe
!   Updated license.
!
!   Revision 1.8  2008/11/19 22:17:34  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.7  2004/02/16 23:49:32  haselbac
!   Changed setting of MOCTR (for PLAG testing) and allocation of NOCTR
!
!   Revision 1.6  2002/11/27 22:34:06  haselbac
!   Changed estimate of MOCTR to resolve occasional error
!
!   Revision 1.5  2002/10/08 15:49:21  haselbac
!   {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
!   Revision 1.4  2002/10/05 19:15:22  haselbac
!   Some cleaning up and fixed bug
!
!   Revision 1.3  2002/09/09 15:11:46  haselbac
!   global now under regions, bug fixes by me and by Tim Baker, pass bounding box
!
!   Revision 1.2  2002/07/25 14:53:03  haselbac
!   Added improved estimate for MOCTR based on Tim Bakers email
!
!   Revision 1.1  2002/06/27 15:48:16  haselbac
!   Initial revision
!
! *****************************************************************************









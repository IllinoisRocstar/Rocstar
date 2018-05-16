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
! Purpose: Modify mesh to conform to boundary displacement using routines 
!   provided by Tim Baker.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   0. Compile this in fixed format, hence use .f extension, otherwise does 
!      not compile...
!   1. Noteworthy changes from Tim Bakers original routines:
!      - Routine SECOND and all calls and associated FORMAT statements 
!        commented out.
!      - In a few places, cleaned out segments of code which Tim commented out.
!   2. Original code assumes implicit typing - AARGH! Changed only routine 
!      RFLU_Repair3D to use IMPLICIT NONE.
!   3. Current routine RFLU_Repair3D was originally PROGRAM Repair3D.
!   4. All routines are restricted to tetrahedra.
!
!******************************************************************************
!
! $Id: RFLU_ModRepair3D.F90,v 1.4 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

MODULE RFLU_ModRepair3D

  INTEGER, PARAMETER :: MNODE  = 100000, &
                        MCELL  = 6*MNODE, & 
                        MEDGE  = 7*MNODE, &
                        MBPTS  = 40000, &
                        MBFACE = 2*MBPTS, &  
                        MOCTR  = 60000,&
                        MTEST  = 50000, & 
                        MCAV   = 24000, &
                        MRING  = 80, &
                        MEDKP  = 2*MRING

  INTEGER :: MXBPTS,MXCAV,MXCELL,MXEDGE,MXFACE,MXNODE,MXOCTR,MXRING,MXTEST

! **************************************************************************** 
! Contained subroutines and functions
! ****************************************************************************
  
      CONTAINS

      SUBROUTINE RFLU_Repair3D(NBPTS,NBFACE,NNODE,NCELL,XI,NFCEI,NDCI, &  
                               XBNDYI,modFlag)

!    ******************************************************************
!    *                                                                *
!    *   THIS PROGRAM MODIFIES A TETRAHEDRAL MESH TO CONFORM WITH A   *
!    *   DISPLACED BOUNDARY SURFACE. THE RESULTING VOLUME MESH WILL,  *
!    *   IN GENERAL, HAVE A MODIFIED TOPOLOGY.                        *
!    *                                                                *
!    ******************************************************************
!
!    ******************************************************************
!    *                                                                *
!    *    COPYRIGHT  (C)  TIM BAKER   2002                            *
!    *                                                                *
!    ******************************************************************
!
!    ******************************************************************
!    *                                                                *
!    *    VERSION 1.1  :  JULY  2002                                  *
!    *                                                                *
!    ******************************************************************
!
!    ******************************************************************


      IMPLICIT NONE

! --- Arguments

      INTEGER, INTENT(INOUT) :: NBFACE,NBPTS,NCELL,NNODE
      INTEGER, INTENT(INOUT) :: NFCEI(3,NBFACE),NDCI(4,NCELL)
      DOUBLE PRECISION, INTENT(INOUT) :: XI(3,NNODE),XBNDYI(3,NBPTS)
      LOGICAL, INTENT(IN) :: modFlag

! --- Locals

      INTEGER :: IOCTR,N,NEDGE,NFAIL
      INTEGER :: IDGP(MNODE),IDONE(MNODE),IEDKP(4,MEDKP),IFLAG(MNODE), &
                 IKEEP(MTEST),IPOINT(MNODE),IPROT(MCELL),IRING(MRING), &
                 ISHK(MRING),HK(MRING), &
                 ITYP(MNODE),KSHK(MRING),KSRCH(MTEST),LDEL(MTEST), &
                 LISTF(MCELL),LNBR(MRING),LNKDN(MCELL),LNKUP(MCELL), & 
                 MNBR(MRING),NACPT(MCELL),NBH(4,MCELL),NBHKP(3,MRING), & 
                 NCAV(4,MCAV),NCAVFC(3,MTEST),NDC(4,MCELL),NDG(2,MEDGE), & 
                 NDGP(MEDGE),NEDGRM(MTEST),NEWC(MRING),NEWCEL(MTEST), &
                 NEWNBH(4,MTEST),NFAD(3,MRING),NFCE(3,MBFACE), &
                 NFILL(MTEST),NFLAG(MCELL),NLINK(MNODE),NOCTR(2,MOCTR), & 
                 NOLD(MTEST),NPOINT(MCELL),NPP(MRING),NPROP(MBFACE), &
                 NPTET(MNODE),NREF(MNODE),NSHAKE(MTEST),NSRCH(MTEST), &
                 NTETKP(MRING),NTRI(3,MTEST),NVCNT(MNODE)
      DOUBLE PRECISION :: RCMX,TOLV,VOLMIN
      DOUBLE PRECISION :: COUNT(MNODE),DENS(MNODE),DS(MRING),DX(MRING), & 
                          DY(MRING),DZ(MRING),FAC(MNODE),RAD(MTEST), &
                          RAT(MCELL),RC(MCELL),RCRIN(MTEST), &
                          RESID(MNODE),SIG1(MCELL),SIG2(MCELL), &
                          SIG3(MCELL),SV(MEDGE),V(MTEST),VLT(MEDKP), &
                          VOL(MCELL),XC(MTEST),XCEN(MCELL),XCH(3,MNODE), & 
                          XFAR(2),XHOLD(2,MTEST),XKEEP(2), &
                          XOCTR(2,MTEST),YC(MTEST),YCEN(MCELL),YFAR(2), &
                          YHOLD(2,MTEST),YKEEP(2),YOCTR(2,MTEST), &
                          ZC(MTEST),ZCEN(MCELL),ZFAR(2),ZHOLD(2,MTEST), &
                          ZOCTR(2,MTEST),ZKEEP(2),X(3,MNODE), &
                          XBNDY(3,MBPTS),XNEWBN(3,MNODE)                         

!    ******************************************************************
!
!    SET TOLERANCES.
!
      TOLV       = 1.E-13

!    SET PARAMETER LIMITS

      MXBPTS     = MBPTS
      MXNODE     = MNODE
      MXCELL     = MCELL
      MXFACE     = MBFACE
      MXEDGE     = MEDGE
      MXOCTR     = MOCTR
      MXCAV      = MCAV
      MXRING     = MRING
      MXTEST     = MTEST

! --- Copy input data into arrays

      DO N = 1,NBFACE
        NFCE(1,N) = NFCEI(1,N)
        NFCE(2,N) = NFCEI(2,N)   
        NFCE(3,N) = NFCEI(3,N)                  
      END DO ! N    

      DO N = 1,NCELL
        NDC(1,N) = NDCI(1,N)
        NDC(2,N) = NDCI(2,N)   
        NDC(3,N) = NDCI(3,N) 
        NDC(4,N) = NDCI(4,N)                          
      END DO ! N    

      DO N = 1,NNODE
        X(1,N) = XI(1,N)
        X(2,N) = XI(2,N)
        X(3,N) = XI(3,N)                
      END DO ! N    
      
      DO N = 1,NBPTS
        XBNDY(1,N) = XBNDYI(1,N)
        XBNDY(2,N) = XBNDYI(2,N)
        XBNDY(3,N) = XBNDYI(3,N)                
      END DO ! N  

!     READ IN SURFACE DATA AND TETRAHEDRAL MESH.
!
!     With the modifications for the incorporation of Repair3D into Rocflu, 
!     the function of INPUT has changed. In the original version of Repair3D, 
!     INPUT read the mesh. In the modified version, it is still called but the
!     mesh is passed through the argument lists. NOTE that ITYP and NPROP are 
!     not passed into RFLU_Repair3D and are still set in INPUT. This will have 
!     to be changed in the future when the boundary triangulation will also
!     change.
!
   5  CALL INPUT (X,NNODE,NDC,NCELL,NFCE,NBPTS,NBFACE, &
                  ITYP,NPROP,XBNDY,XFAR,YFAR,ZFAR)
!
!
!
      CALL STRUCT (X,NNODE,NDC,NBH,IPROT,NCELL,NFCE,NBFACE, &
                   NEDGE,NDG,IDGP,NDGP,IPOINT,NPOINT,NPTET, &
                   XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                   XFAR,YFAR,ZFAR, &
                   IOCTR,NLINK,NOCTR,IDONE,NREF,VOLMIN,RCMX)
!
!     COMPUTE STATISTICS OF INITIAL MESH.
!
      CALL RADRAT (X,NNODE,NDC,NBH,IPROT,NCELL,NFCE,NBFACE, &
                   ITYP,IPOINT,VOL,RC,RAT)
!
!     DISPLACE MESH TO CONFORM WITH NEW BOUNDARY POSITION.
!
      CALL TETMV (X,XNEWBN,NNODE,NDC,NBH,NCELL,NFCE,NBFACE,NFAIL, &
                  ITYP,XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                  NEDGE,NDG,IDGP,NDGP,IPOINT,COUNT,XCH,RESID,SV, &
                  SIG1,SIG2,SIG3,XFAR,YFAR,ZFAR,TOLV)

      IF (NFAIL.NE.0) THEN
        GO TO 100
      ENDIF
!
!     MODIFY MESH TOPOLOGY TO IMPROVE QUALITY OF TETRAHEDRAL
!
      IF ( modFlag .EQV. .TRUE. ) THEN 
        CALL TETMOD (X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                     NDG,IDGP,NDGP,NEDGE,NFCE,NBFACE, &
                     VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                     SIG1,SIG2,SIG3,NVCNT,RESID,COUNT,FAC, &
                     IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                     XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                     XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                     IPOINT,NPOINT,IFLAG,NFLAG, &
                     DX,DY,DZ,DS,VLT,IRING,NTETKP,NFAD,NEWC, &
                     NBHKP,IEDKP,LNBR,ISHK,MNBR,KSHK,NPP, &
                     NFILL,NEWCEL,NTRI, &
                     NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                     NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                     LISTF,VOLMIN,RCMX,TOLV)
!
!       COMPUTE STATISTICS OF NEW MESH.
!
        CALL RADRAT (X,NNODE,NDC,NBH,IPROT,NCELL,NFCE,NBFACE, &
                     ITYP,IPOINT,VOL,RC,RAT)
      END IF ! modFlag

      STOP
!
!     ALGORITHM HAS FAILED IN MESH MOVEMENT PHASE. WRITE OUT
!     MESH DATA FOR THE CORRUPTED MESH.
!
  100 WRITE (6,610)
      STOP
  610 FORMAT(//5X,'FAILURE IN MESH MOVEMENT PHASE')
  620 FORMAT(/5X,'****************************************************'/ &
              5X,'****************************************************'/ &
              5X,'***                                              ***'/ &
              5X,'** CYCLE ',I3,' OF MESH MODIFICATION IS COMPLETE  **'/ &
              5X,'**          PARAMETRIC TIME = ',F8.4,'            **'/ &
              5X,'***                                              ***'/ &
              5X,'****************************************************'/ &
              5X,'****************************************************')
  630 FORMAT(/5X,'MESH MOVEMENT SCHEME HAS FAILED IN STEP ',I4,' OF A ', &
              I4,' STAGE CYCLE.'/ &
              5X,'DOUBLE THE NUMBER OF STAGES AND REPEAT PROCEDURE.'/) 
     
      END SUBROUTINE RFLU_Repair3D
      

      

!
!     ******************************************************************
!
      SUBROUTINE INPUT (X,NNODE,NDC,NCELL,NFCE,NBPTS,NBFACE, &
                        ITYP,NPROP,XBNDY,XFAR,YFAR,ZFAR)
!
!     ******************************************************************
!     *                                                                *
!     *    READ IN SURFACE AND VOLUME MESHES                           *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1998                            *
!     *                                                                *
!     ******************************************************************
!
!
      IMPLICIT NONE

      INTEGER :: NBFACE,NBPTS,NCELL,NNODE
      INTEGER :: NDC(4,*)
      INTEGER :: ITYP(*),NPROP(*)
      INTEGER :: NFCE(3,*)
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: XBNDY(3,*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)
      
      INTEGER :: J,JAMIN,JMAX,JMIN,L,N,N1,N2,N3
      DOUBLE PRECISION :: ANGL,ANGL1,ANGL2,ANGL3,ANGMAX,ANGMIN,FAREA, &
                          FMAX,FMIN,Q,QMAX,QMIN
!
!     ******************************************************************
!
!      IREAD      = 5
!
!     READ IN TETRAHEDRAL MESH
!
!      READ (IREAD,500)
!      READ (IREAD,520) NBPTS,NBFACE,NNODE,NCELL
      
! --- Check sizes of arrays against parameters      
      
      IF ( NBPTS > MXBPTS ) THEN 
        GO TO 300
      END IF ! NBPTS
      
      IF ( NNODE > MXNODE ) THEN 
        GO TO 305
      END IF ! NNODE
      
      IF ( NBFACE > MXFACE ) THEN 
        GO TO 310
      END IF ! NBFACE
      
      IF ( NCELL > MXCELL ) THEN 
        GO TO 315
      END IF ! NCELL
      
! --- Read coordinates, boundary faces, tetrahedra connectivity      
      
!      READ (IREAD,500)
!      DO 10 N=1,NNODE
!      READ (IREAD,515) X(1,N),X(2,N),X(3,N),ITYP(N)
!   10 CONTINUE
!      READ (IREAD,500)
!      DO 15 L=1,NBFACE
!      READ (IREAD,520) (NFCE(K,L),K=1,3),NPROP(L)
!   15 CONTINUE
!      READ (IREAD,500)
!      DO 25 L=1,NCELL
!      READ (IREAD,520) (NDC(K,L),K=1,4)
!   25 CONTINUE

! --- ITYP and NPROP must be initialized because they are not passed 
!     into RFLU_Repair3D

      DO N = 1,NNODE
        ITYP(N) = 1
      END DO ! N

      DO L = 1,NBFACE
        NPROP(L) = 1
      END DO ! L

      WRITE (6,700) NBFACE,NNODE,NCELL

!     COMPUTE AREA AND ANGLES OF EACH SURFACE TRIANGLE

      JMAX        = 1
      JMIN        = 1
      JAMIN       = 1
      ANGMIN      = 360.
      ANGMAX      = 0.
      QMIN        = 1.E15
      QMAX        = 1.
      DO 75 J=1,NBFACE
      N1          = NFCE(1,J)
      N2          = NFCE(2,J)
      N3          = NFCE(3,J)

      CALL FANGLE (J,X,NFCE,ANGL1,ANGL2,ANGL3,Q)

      ANGL        = MIN(ANGL1,ANGL2,ANGL3)
      IF (ANGL.GT.ANGMIN) GO TO 60
      JAMIN       = J
      ANGMIN      = ANGL
   60 ANGMAX      = MAX(ANGL1,ANGL2,ANGL3,ANGMAX)
      QMIN        = MIN(Q,QMIN)
      QMAX        = MAX(Q,QMAX)
      FAREA       = FACEAR (X,N1,N2,N3)
      IF (J.EQ.1) FMAX = FAREA
      IF (J.EQ.1) FMIN = FAREA
      IF (FAREA.LT.FMAX) GO TO 70
      JMAX        = J
      FMAX        = FAREA
   70 IF (FAREA.GT.FMIN) GO TO 75
      JMIN        = J
      FMIN        = FAREA
   75 CONTINUE
      WRITE (6,705) ANGMIN,ANGMAX,QMIN,QMAX
      N1          = NFCE(1,JAMIN)
      N2          = NFCE(2,JAMIN)
      N3          = NFCE(3,JAMIN)
      IF (ANGMIN.LT.1.) WRITE (6,706) JAMIN,N1,X(1,N1),X(2,N1),X(3,N1), &
                                            N2,X(1,N2),X(2,N2),X(3,N2), &
                                            N3,X(1,N3),X(2,N3),X(3,N3)
      N1          = NFCE(1,JMIN)
      N2          = NFCE(2,JMIN)
      N3          = NFCE(3,JMIN)
      WRITE (6,707) FMAX,FMIN,JMIN,N1,X(1,N1),X(2,N1),X(3,N1), &
                                   N2,X(1,N2),X(2,N2),X(3,N2), &
                                   N3,X(1,N3),X(2,N3),X(3,N3)

!     DETERMINE MINIMUM AND MAXIMUM EXTENT OF GEOMETRY

      XFAR(1)    = X(1,1)
      XFAR(2)    = X(1,1)
      YFAR(1)    = X(2,1)
      YFAR(2)    = X(2,1)
      ZFAR(1)    = X(3,1)
      ZFAR(2)    = X(3,1)
      DO 90 N=2,NNODE
      XFAR(1)    = MIN(XFAR(1),X(1,N))
      XFAR(2)    = MAX(XFAR(2),X(1,N))
      YFAR(1)    = MIN(YFAR(1),X(2,N))
      YFAR(2)    = MAX(YFAR(2),X(2,N))
      ZFAR(1)    = MIN(ZFAR(1),X(3,N))
      ZFAR(2)    = MAX(ZFAR(2),X(3,N))
   90 CONTINUE
      WRITE (6,920) XFAR(1),XFAR(2),YFAR(1),YFAR(2),ZFAR(1),ZFAR(2)
  920 FORMAT('XFAR ',F10.2,1X,F10.2,' YFAR ',F10.2,1X,F10.2, &
             ' ZFAR ',F10.2,1X,F10.2)
!
!     READ IN COORDINATES OF NEW SURFACE POINT POSITIONS
!
!      IDATA     = 10
!      READ (IDATA,500)
!      READ (IDATA,510) NBPTS
!      READ (IDATA,500)
!      DO 100 N=1,NBPTS
!      READ (IDATA,525) XBNDY(1,N),XBNDY(2,N),XBNDY(3,N)
!  100 CONTINUE
  
!      TIM        = SECOND (0)
!      WRITE (6,710) TIM

      RETURN
  300 WRITE (6,600) NBPTS
      STOP
  305 WRITE (6,605) NNODE
      STOP
  310 WRITE (6,610) NBFACE
      STOP
  315 WRITE (6,615) NCELL
      STOP
  500 FORMAT(1X)
  510 FORMAT(I10)
! 515 FORMAT(3F12.5,I10)
! 515 FORMAT(3F13.5,I10)
! 515 FORMAT(3F10.4,I10)
  515 FORMAT(3E13.5,I10)
  520 FORMAT(6I10)
  525 FORMAT(3E13.5)
  600 FORMAT(///5X,'NUMBER OF SURFACE POINTS ',I7/ &
                5X,'EXCEEDS MAXIMUM ALLOWED. INCREASE SIZE OF MBPTS'/ &
                5X,'PROGRAM STOPPED IN ROUTINE INPUT')
  605 FORMAT(///5X,'TOTAL NUMBER OF MESH POINTS ',I7/ &
                5X,'EXCEEDS MAXIMUM ALLOWED. INCREASE SIZE OF MNODE'/ &
                5X,'PROGRAM STOPPED IN ROUTINE INPUT')
  610 FORMAT(///5X,'NUMBER OF SURFACE FACES ',I7/ &
                5X,'EXCEEDS MAXIMUM ALLOWED. INCREASE SIZE OF MBPTS'/ &
                5X,'PROGRAM STOPPED IN ROUTINE INPUT')
  615 FORMAT(///5X,'TOTAL NUMBER OF MESH CELLS ',I7/ &
                5X,'EXCEEDS MAXIMUM ALLOWED. INCREASE SIZE OF MCELL'/ &
                5X,'PROGRAM STOPPED IN ROUTINE INPUT')
  700 FORMAT(/5X,'SURFACE AND VOLUME MESH READ'// &
              5X,'NUMBER OF BOUNDARY SURFACE TRIANGLES      = ',I7/ &
              5X,'TOTAL NUMBER OF MESH POINTS               = ',I7/ &
              5X,'TOTAL NUMBER OF MESH CELLS                = ',I7)
  705 FORMAT(/5X,'MINIMUM BOUNDARY FACE ANGLE = ',F6.2/ &
              5X,'MAXIMUM BOUNDARY FACE ANGLE = ',F6.2/ &
              5X,'MINIMUM BOUNDARY RADIUS RATIO = ',F6.2/ &
              5X,'MAXIMUM BOUNDARY RADIUS RATIO = ',F6.2)
  706 FORMAT(/5X,'WARNING !   MINIMUM FACE ANGLE IS LESS THAN 1 DEGREE'/ &
              5X,'FACE ADDRESS ',I6/ &
              5X,'VERTEX ',I6,' X = ',F12.4,' Y = ',F12.4,' Z = ',F12.4/ &
              5X,'VERTEX ',I6,' X = ',F12.4,' Y = ',F12.4,' Z = ',F12.4/ &
              5X,'VERTEX ',I6,' X = ',F12.4,' Y = ',F12.4,' Z = ',F12.4)
  707 FORMAT(/5X,'MAX BOUNDARY FACE AREA = ',E13.5/ &
              5X,'MIN BOUNDARY FACE AREA = ',E13.5/ &
              5X,'ADDRESS OF BOUNDARY FACE WITH MINIMUM AREA',I6/ &
              5X,'VERTEX ',I6,' X = ',F12.4,' Y = ',F12.4,' Z = ',F12.4/ &
              5X,'VERTEX ',I6,' X = ',F12.4,' Y = ',F12.4,' Z = ',F12.4/ &
              5X,'VERTEX ',I6,' X = ',F12.4,' Y = ',F12.4,' Z = ',F12.4/)
      END SUBROUTINE INPUT





!
!     ******************************************************************
!
      SUBROUTINE STRUCT (X,NNODE,NDC,NBH,IPROT,NCELL,NFCE,NBFACE, &
                         NEDGE,NDG,IDGP,NDGP,IPOINT,NPOINT,NPTET, &
                         XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                         XFAR,YFAR,ZFAR, &
                         IOCTR,NLINK,NOCTR,IDONE,NREF,VOLMIN,RCMX) 
!
!     ******************************************************************
!     *                                                                *
!     *  CREATE DATA STRUCTURE FOR USE IN VOLUME MESH                  *
!     *  RECONSTRUCTION  BASED ON DELAUNAY REFINEMENT                  *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1998                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,NBFACE,NCELL,NEDGE,NNODE
      INTEGER :: NDC(4,*),NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: NBH(4,*),IPROT(*),NPTET(*),NPOINT(*),IPOINT(*)
      INTEGER :: NFCE(3,*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)      
      DOUBLE PRECISION :: VOLMIN,RCMX
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)

      INTEGER :: I,IEDG,J,K,KK,L,LK,LM,LN,L1,L2,MM,M1,M2,M3,M4,N,NCELM, &
                 NCHG,NCNT,NEND,NINC,NLOOK,NN,NOLD,N1,N2,N3,N4
      DOUBLE PRECISION :: AREA,C5,C6,C9,FAC,H1,H2,H3,H4,RNX,RNY,RNZ,TOL, &
                          VCELL,VMAX,VMIN,VRAT,VSEP,XSHF,YSHF,ZSHF
      DOUBLE PRECISION :: V(4),AR(4)
!
!     ******************************************************************
!
      TOL       = 1.0E-9
!
!     INITIALIZE NPTET AND IPOINT ARRAY
!
      DO 10 N=1,NNODE
      IPOINT(N) = 0
      NPTET(N)  = 0
   10 CONTINUE
!
!     COMPUTE VOLUME, CIRCUMRADIUS AND RADIUS RATIO
!
      DO 40 L=1,NCELL
      N1        = NDC(1,L)
      N2        = NDC(2,L)
      N3        = NDC(3,L)
      N4        = NDC(4,L)
      NEND      = N4
      NCNT      = 0
   15 NCNT      = NCNT  +1
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
!
!     COMPUTE VOLUME OF TETRAHEDRON AND AREAS ALL FOUR FACES
!
      V(NCNT) = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                         +RNZ*(X(3,N4)  -X(3,N1))
      AR(NCNT) = SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      IF (N1.EQ.NEND) GO TO 20
      N       = N1
      N1      = N2
      N2      = N3
      N3      = N4
      N4      = N
      GO TO 15
   20 VOL(L)  = .25*(ABS(V(1))  +ABS(V(2))  +ABS(V(3))  +ABS(V(4)))
      AREA    = AR(1)  +AR(2)  +AR(3)  +AR(4)
!
!     COMPUTE CIRCUMCENTER COORDINATES
!
      RC(L)     = 0.
      RAT(L)    = 1.E10
      IF (VOL(L).LT.TOL) GO TO 30
      FAC       = .5/V(4)
      XSHF      = .5*(X(1,N1)  +X(1,N2)  +X(1,N3)  +X(1,N4))
      YSHF      = .5*(X(2,N1)  +X(2,N2)  +X(2,N3)  +X(2,N4))
      ZSHF      = .5*(X(3,N1)  +X(3,N2)  +X(3,N3)  +X(3,N4))
      H1        = (X(1,N4)  -X(1,N1))*(X(1,N4)  -XSHF  +X(1,N1)) &
                 +(X(2,N4)  -X(2,N1))*(X(2,N4)  -YSHF  +X(2,N1)) &
                 +(X(3,N4)  -X(3,N1))*(X(3,N4)  -ZSHF  +X(3,N1))
      H2        = (X(1,N3)  -X(1,N1))*(X(1,N3)  -XSHF  +X(1,N1)) &
                 +(X(2,N3)  -X(2,N1))*(X(2,N3)  -YSHF  +X(2,N1)) &
                 +(X(3,N3)  -X(3,N1))*(X(3,N3)  -ZSHF  +X(3,N1))
      H3        = (X(1,N2)  -X(1,N1))*(X(1,N2)  -XSHF  +X(1,N1)) &
                 +(X(2,N2)  -X(2,N1))*(X(2,N2)  -YSHF  +X(2,N1)) &
                 +(X(3,N2)  -X(3,N1))*(X(3,N2)  -ZSHF  +X(3,N1))
      C5        = H2*(X(3,N2)  -X(3,N1))  -H3*(X(3,N3)  -X(3,N1))
      C6        = H2*(X(2,N2)  -X(2,N1))  -H3*(X(2,N3)  -X(2,N1))
      C9        = H3*(X(1,N3)  -X(1,N1))  -H2*(X(1,N2)  -X(1,N1))
      XCEN(L)   = .5*XSHF  +(H1*RNX  +(X(2,N4)  -X(2,N1))*C5 &
                                     -(X(3,N4)  -X(3,N1))*C6)*FAC
      YCEN(L)   = .5*YSHF  +(-(X(1,N4)  -X(1,N1))*C5  +H1*RNY &
                             -(X(3,N4)  -X(3,N1))*C9)*FAC
      ZCEN(L)   = .5*ZSHF  +((X(1,N4)  -X(1,N1))*C6 &
                            +(X(2,N4)  -X(2,N1))*C9  +H1*RNZ)*FAC
      RC(L)     = SQRT((X(1,N1)  -XCEN(L))**2  +(X(2,N1)  -YCEN(L))**2 &
                                               +(X(3,N1)  -ZCEN(L))**2)
      RAT(L)    = RC(L)*AREA/VOL(L)
      VMIN      = MIN(ABS(V(1)),ABS(V(2)),ABS(V(3)),ABS(V(4)))
      VMAX      = MAX(ABS(V(1)),ABS(V(2)),ABS(V(3)),ABS(V(4)))
      VSEP      = VMAX  -VMIN
      VRAT      = VSEP/VOL(L)
!     IF (VSEP.GT.TOL) WRITE (6,610) VSEP,VOL(L),VRAT
      IF (L.EQ.1) VOLMIN = VOL(1)
      IF (L.EQ.1) RCMX = RC(1)
      VOLMIN    = MIN(VOLMIN,VOL(L))
      RCMX      = MAX(RCMX,RC(L))
      GO TO 40
   30 VCELL     = VOL(L)/6.
      WRITE (6,600) VCELL
   40 CONTINUE
!
!     SET UP GHOST CELLS BEHIND BOUNDARY FACES
!
      NCELM     = NCELL
      DO 50 J=1,NBFACE
      NCELL     = NCELL  +1
      NDC(1,NCELL) = NFCE(1,J)
      NDC(2,NCELL) = NFCE(2,J)
      NDC(3,NCELL) = NFCE(3,J)
      NDC(4,NCELL) = -1
   50 CONTINUE
!
!     INITIALIZE NPOINT AND SET STARTING VALUES FOR NPTET AND IPROT
!
      DO 55 L=1,NCELL
      NPOINT(L) = 0
      IPROT(L)  = L
      DO 55 K=1,4
      N         = NDC(K,L)
      IF (N.LT.0) GO TO 55
      IF (NPTET(N).GT.0) GO TO 55
      NPTET(N)  = L
   55 CONTINUE
!
!     CREATE NBH ARRAY CONTAINING ADDRESSES OF NEIGHBORING TETRAHEDRA
!
      NLOOK     = NCELL
   60 DO 120 LK=1,NLOOK
      L         = IPROT(LK)
      IF (NPOINT(L).EQ.4) GO TO 120
      N1        = NDC(1,L)
      N2        = NDC(2,L)
      N3        = NDC(3,L)
      N4        = NDC(4,L)
      IF (N1.LT.0) GO TO 100
   65 J         = NPTET(N1)
      IF (J.EQ.0.OR.J.EQ.L) GO TO 100
      M1        = NDC(1,J)
      M2        = NDC(2,J)
      M3        = NDC(3,J)
      M4        = NDC(4,J)
   70 IF (M1.EQ.N1) GO TO 75
      MM        = M1
      M1        = M2
      M2        = M3
      M3        = M4
      M4        = MM
      GO TO 70
   75 NINC      = 0
      IF (N2.EQ.M2.OR.N2.EQ.M3.OR.N2.EQ.M4) NINC = NINC  +1
      IF (N3.EQ.M2.OR.N3.EQ.M3.OR.N3.EQ.M4) NINC = NINC  +1
      IF (N4.EQ.M2.OR.N4.EQ.M3.OR.N4.EQ.M4) NINC = NINC  +1
      IF (NINC.LT.2) GO TO 100
   80 IF (NPOINT(L).EQ.0) GO TO 90
      DO 85 K=1,NPOINT(L)
      LN        = NBH(K,L)
      IF (LN.EQ.J) GO TO 100
   85 CONTINUE
   90 NPOINT(L) = NPOINT(L)  +1
      NBH(NPOINT(L),L) = J
      NPOINT(J) = NPOINT(J)  +1
      NBH(NPOINT(J),J) = L
  100 IF (N1.EQ.NDC(4,L)) GO TO 120
      NN        = N1
      N1        = N2
      N2        = N3
      N3        = N4
      N4        = NN
      IF (N1.LT.0) GO TO 100
      GO TO 65
  120 CONTINUE
!
!     RESET IPROT AND NPTET PRIOR TO REPEAT SEARCH FOR CELL NEIGHBORS
!
      NCHG      = 0
      DO 125 LK=1,NLOOK
      L         = IPROT(LK)
      IF (NPOINT(L).EQ.4) GO TO 125
      NCHG      = NCHG  +1
      IPROT(NCHG) = L
  125 CONTINUE
      IF (NCHG.EQ.0) GO TO 140
      NLOOK     = NCHG
      DO 130 LK=1,NLOOK
      L         = IPROT(LK)
      DO 130 K=1,4
      N         = NDC(K,L)
      IF (N.GT.0) NPTET(N) = 0
  130 CONTINUE
      DO 135 LK=1,NLOOK
      L         = IPROT(LK)
      DO 135 K=1,4
      N         = NDC(K,L)
      IF (N.LT.0) GO TO 135
      IF (NPTET(N).EQ.0) NPTET(N)  = L
  135 CONTINUE
      GO TO 60
!
!     ALL CELL NEIGHBORS HAVE BEEN FOUND. CHECK WHETHER EVERY CELL HAS
!     EXACTLY FOUR NEIGHBORS. THEN RE-INITIALIZE NPOINT ARRAY AND SET
!     VALUES OF IPROT FOR USE BY THE MESH RECONSTRUCTION ALGORITHM
!
  140 DO 155 L=1,NCELL
      IF (NPOINT(L).NE.4) GO TO 320
      NPOINT(L) = 0
      IPROT(L)  = 1
      IF (L.GT.NCELM) GO TO 155
      IPROT(L)  = 0
      DO 150 K=1,4
      N         = NDC(K,L)
      IF (N.GT.0) NPTET(N) = L
  150 CONTINUE
  155 CONTINUE
!
!     INITIALIZE OCTREE STRUCTURE
!
      IOCTR      = 1
      NLINK(1)   = 0
      NOCTR(1,1) = 1
      NOCTR(2,1) = 0
      IDONE(1)   = 1
      NREF(1)    = 1
      IPOINT(1)  = 0
      DO 170 N=2,NNODE

      CALL OCTFIL (N,X,NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR)

      IDONE(N)   = 1
      IPOINT(N)  = 0
  170 CONTINUE
!
!     CREATE EDGE DATA STRUCTURE FOR MESH
!
      DO 175 N=1,NNODE
      IDGP(N)    = 0
  175 CONTINUE
      NEDGE      = 0
      DO 200 J=1,NCELL
      IF (NDC(4,J).EQ.-1) GO TO 200
      I          = 1
  180 I          = I  +1
      IF (I.EQ.5) GO TO 200
      N1         = NDC(I,J)
      K          = 0
  185 K          = K  +1
      IF (K.EQ.I) GO TO 180
      N2         = NDC(K,J)
      L1         = MIN(N1,N2)
      L2         = MAX(N1,N2)
      NOLD       = IDGP(L1)
  190 IF (NOLD.EQ.0) GO TO 195
      IF (L2.EQ.MAX(NDG(1,NOLD),NDG(2,NOLD))) GO TO 185
      IEDG       = NOLD
      NOLD       = NDGP(NOLD)
      GO TO 190
  195 NEDGE      = NEDGE  +1
      NDG(1,NEDGE) = L1
      NDG(2,NEDGE) = L2
      NDGP(NEDGE)  = 0
      IF (IDGP(L1).NE.0) NDGP(IEDG) = NEDGE
      IF (IDGP(L1).EQ.0) IDGP(L1)   = NEDGE
      GO TO 185
  200 CONTINUE
      RETURN
  320 WRITE (6,620) L,(NDC(KK,L),KK=1,4),NPOINT(L)
      STOP
  600 FORMAT(//5X,'TETRAHEDRON WITH AN EXTREMELY SMALL VOLUME FOUND'// &
                  ' IN ROUTINE STRUCT'//5X,'VOLUME = ',E13.5/)
  610 FORMAT(5X,'IMPRECISE ESTIMATE OF TETRAHEDRON VOLUME', &
             5X,'VSEP = ',E13.5,' VOLUME = ',E13.5, &
                ' VSEP/VOLUME = ',E13.5)
  620 FORMAT(//5X,'NPOINT IS NOT EQUAL TO FOUR FOR TETRAHEDRON ',I7, &
               5X,'TETRAHEDRON VERTICES ARE ',4I6,' AND NPOINT = ',I2, &
               5X,'PROGRAM STOPPED IN ROUTINE STRUCT')
      END SUBROUTINE STRUCT



!
!     ******************************************************************
!
      SUBROUTINE TETMV (X,XNEWBN,NNODE,NDC,NBH,NCELL,NFCE,NBFACE,NFAIL, &
                        ITYP,XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                        NEDGE,NDG,IDGP,NDGP,IPOINT,FCOUNT,XC,RESID,SV, &
                        SIG1,SIG2,SIG3,XFAR,YFAR,ZFAR,TOLV) 

!     ******************************************************************
!     *                                                                *
!     *   LAPLACIAN SOLVER TO RELAX A TETRAHEDRAL MESH TO CONFORM      *
!     *   WITH A NEW BOUNDARY SURFACE POSITION.                        *
!     *   THE LAPLACIAN APPROXIMATION IS BASED ON THE ACCUMMULATION    *
!     *   OF EDGE DIFFERENCES AND WILL IN GENERAL FAIL TO GIVE ZERO    *
!     *   FOR A LINEAR VARIATION IN DISPLACEMENT. THIS ROUTINE WILL    *
!     *   THEREFORE ONLY BE EFFECTIVE FOR SMALL PERTURBATIONS OF THE   *
!     *   BOUNDARY SURFACE.                                            *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1999                            *
!     *                                                                *
!     ******************************************************************


!     ******************************************************************
!     *                                                                *
!     *    DIMENSION CONTROL PARAMETERS                                *
!     *                                                                *
!     ******************************************************************

      IMPLICIT NONE

      INTEGER :: NBFACE,NCELL,NEDGE,NFAIL,NNODE
      INTEGER :: NDC(4,*),NBH(4,*),ITYP(*),NFCE(3,*)
      INTEGER :: NDG(2,*),IDGP(*),NDGP(*),IPOINT(*)      
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*),XNEWBN(3,*)
      DOUBLE PRECISION :: XCEN(*),YCEN(*),ZCEN(*),VOL(*),RC(*),RAT(*)
      DOUBLE PRECISION :: FCOUNT(*),XC(3,*),RESID(*),SV(*), &
                          SIG1(*),SIG2(*),SIG3(*)
      DOUBLE PRECISION :: RMAX(3),ITMAX(3)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)
      DOUBLE PRECISION :: XCN,YCN,ZCN,VL,RADC
      DOUBLE PRECISION :: VMIN,VMAX
      
      INTEGER :: I,IBIG,IEDG,ISMALL,IT,J,K,KP1,KP2,KP3,L,LOOP,LSUB, &
                 L1,L2,N,NNEG,NOLD,NPASS,N1,N2,N3,N4
      DOUBLE PRECISION :: AREA,CNDMAX,CNDMIN,DIFX,EPS,FAC,RBIG,RESM, &
                          RMAX1,RNX,RNY,RNZ,SIGMX,SIGMN,STRMAX,STRMIN, &
                          STRTET,THIRD,VTET,VTET1,VTET2,XBAR,YBAR,ZBAR
!
!     ******************************************************************
!
      NFAIL     = 0
      THIRD     = 1./3.
      EPS       = 1.
!
!     SUB-ITERATE TO MOVE MESH IN LSUB STEPS
!
      LSUB      = 1
      FAC       = 1./FLOAT(LSUB)
      DO 150 L=1,LSUB
!
!     INITIALIZE VERTEX EDGE COUNT, IPOINT AND IDGP ARRAY
!
      DO 10 N=1,NNODE
      FCOUNT(N)  = 0.
      IPOINT(N)  = 0
      IDGP(N)    = 0
   10 CONTINUE
!
!     COMPUTE VOLUME OF EACH TETRAHEDRON
!
      DO 22 J=1,NCELL
      IF (NBH(1,J).EQ.0) GO TO 22
      IF (NDC(4,J).EQ.-1) GO TO 22
      VOL(J)     = 0.
      DO 20 K=1,4
      KP1        = MOD(K,4)  +1
      KP2        = MOD(KP1,4)  +1
      KP3        = MOD(KP2,4)  +1
      N1         = NDC(K,J)
      N2         = NDC(KP1,J)
      N3         = NDC(KP2,J)
      N4         = NDC(KP3,J)
      XBAR      = THIRD*(X(1,N1)  +X(1,N2)  +X(1,N3))
      YBAR      = THIRD*(X(2,N1)  +X(2,N2)  +X(2,N3))
      ZBAR      = THIRD*(X(3,N1)  +X(3,N2)  +X(3,N3))
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      VTET      = (X(1,N4)  -XBAR)*RNX  +(X(2,N4)  -YBAR)*RNY &
                 +(X(3,N4)  -ZBAR)*RNZ
      VOL(J)    = VOL(J)  +ABS(VTET)
   20 CONTINUE
   22 CONTINUE
      NPASS     = 0
      DO 25 J=1,NCELL
      IF (NBH(1,J).EQ.0) GO TO 25
      IF (NDC(4,J).EQ.-1) GO TO 25
      IF (NPASS.EQ.0) THEN
         NPASS     = 1
         VMIN      = VOL(J)
         VMAX      = VOL(J)
      ELSE
         VMIN      = MIN(VMIN,VOL(J))
         VMAX      = MAX(VMAX,VOL(J))
      ENDIF
   25 CONTINUE
      DO 30 J=1,NCELL
      IF (NBH(1,J).EQ.0) GO TO 30
      IF (NDC(4,J).EQ.-1) GO TO 30
      VOL(J)    = 1.  +(VMAX  -VMIN)/VOL(J)
   30 CONTINUE
      WRITE (6,940) VMIN,VMAX
  940 FORMAT(/'VMIN = ',E13.5,'   VMAX = ',E13.5/)
!
!     CREATE EDGE DATA STRUCTURE FOR MESH
!
      NEDGE      = 0
      DO 70 J=1,NCELL
      IF (NBH(1,J).EQ.0) GO TO 70
      IF (NDC(4,J).EQ.-1) GO TO 70
      I          = 1
   50 I          = I  +1
      IF (I.EQ.5) GO TO 70
      N1         = NDC(I,J)
      K          = 0
   55 K          = K  +1
      IF (K.EQ.I) GO TO 50
      N2         = NDC(K,J)
      L1         = MIN(N1,N2)
      L2         = MAX(N1,N2)
      NOLD       = IDGP(L1)
   60 IF (NOLD.EQ.0) GO TO 65
      IF (L2.EQ.MAX(NDG(1,NOLD),NDG(2,NOLD))) THEN
         SV(NOLD)    = SV(NOLD)  +VOL(J)
         GO TO 55
      ENDIF
      IEDG       = NOLD
      NOLD       = NDGP(NOLD)
      GO TO 60
   65 NEDGE      = NEDGE  +1
      NDG(1,NEDGE) = L1
      NDG(2,NEDGE) = L2
      NDGP(NEDGE)  = 0
      IF (IDGP(L1).NE.0) NDGP(IEDG) = NEDGE
      IF (IDGP(L1).EQ.0) IDGP(L1)   = NEDGE
      SV(NEDGE)  = VOL(J)
      GO TO 55
   70 CONTINUE
!
!     LOOP OVER X, Y AND Z COMPONENTS OF MESH DISPLACEMENT
!
      LOOP       = 0
   80 LOOP       = LOOP  +1
      IT         = 0
!
!     INITIALIZE POINT DISPLACEMENTS
!
      DO 85 N=1,NNODE
      IPOINT(N)  = 0
      IF (ITYP(N).LT.0) GO TO 85
      XC(LOOP,N) = 0.
   85 CONTINUE
!
!     COMPUTE BOUNDARY POINT DISPLACEMENTS
!
      DO 90 I=1,NBFACE
      DO 90 K=1,3
      N          = NFCE(K,I)
      IF (N.LT.0) GO TO 90
      IF (ITYP(N).LT.0) GO TO 90
      IF (IPOINT(N).GT.0) GO TO 90
      IPOINT(N)  = 1
      XC(LOOP,N) = (XNEWBN(LOOP,N)  -X(LOOP,N))*L*FAC
   90 CONTINUE
      IF (LOOP.GT.1) GO TO 105
!
!     COMPUTE VERTEX DEGREE AND ITS RECIPROCAL AT EACH NODE
!
      DO 95 I=1,NEDGE
      N1         = NDG(1,I)
      N2         = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 95
      IF (ITYP(N1).LT.0.OR.ITYP(N2).LT.0) GO TO 95
      FCOUNT(N1) = FCOUNT(N1)  +SV(I)
      FCOUNT(N2) = FCOUNT(N2)  +SV(I)
   95 CONTINUE
      DO 100 N=1,NNODE
      IF (ITYP(N).LT.0) GO TO 100
      FCOUNT(N)  = 1./FCOUNT(N)
  100 CONTINUE
!
!     START OF ITERATIVE CYCLE
!
  105 IT         = IT  +1
!
!     SET INITIAL RESIDUALS TO ZERO
!
      DO 110 N=1,NNODE
      RESID(N)   = 0.
  110 CONTINUE
!
!     ACCUMMULATE EDGE DIFFERENCES OF DISPLACEMENT COMPONENT FOR ALL
!     EDGES INCIDENT TO EACH POINT
!
      DO 115 I=1,NEDGE
      N1         = NDG(1,I)
      N2         = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 115
      IF (ITYP(N1).LT.0.OR.ITYP(N2).LT.0) GO TO 115
      DIFX       = (XC(LOOP,N1)  -XC(LOOP,N2))*SV(I)
      RESID(N1)  = RESID(N1)  -DIFX
      RESID(N2)  = RESID(N2)  +DIFX
  115 CONTINUE
!
!     UPDATE DISPLACEMENT COMPONENT AT EACH NON-BOUNDARY POINT 
!
      RMAX(LOOP) = 0.
      DO 120 N=1,NNODE
      IF (ITYP(N).LT.0) GO TO 120
      IF (IPOINT(N).EQ.0) THEN
         XC(LOOP,N) = XC(LOOP,N)  +EPS*FCOUNT(N)*RESID(N)
         RESM       = ABS(RESID(N))
         RMAX(LOOP) = MAX(RMAX(LOOP),FCOUNT(N)*RESM)
      ENDIF
  120 CONTINUE

!     IF (IT.EQ.1.OR.MOD(IT,40).EQ.0) WRITE (6,900) IT,RMAX(LOOP)
! 900 FORMAT('IN TETMV...,  IT = ',I4,' RMAX = ',F10.6)

      IF (IT.EQ.1) RMAX1 = RMAX(LOOP)
      IF (RMAX(LOOP).GT.MAX(0.0001*RMAX1,1.0D0-9).AND.IT.LT.2000) &
        GO TO 105
      ITMAX(LOOP) = IT
      IF (LOOP.LT.3) GO TO 80
      RBIG       = MAX(RMAX(1),RMAX(2),RMAX(3))
      IBIG       = MAX(ITMAX(1),ITMAX(2),ITMAX(3))
      WRITE (6,600) RBIG,IBIG
!
!     ITERATIVE CYCLE IS COMPLETE. CHECK FOR SIGN CHANGE IN THE VOLUME
!     OF ANY TETRAHEDRON
!
      STRMAX  = 0.
      STRMIN  = 1.E6
      NNEG    = 0
      DO 130 J=1,NCELL
      IF (NBH(1,J).EQ.0) GO TO 130
      IF (NDC(4,J).EQ.-1) GO TO 130
      N1      = NDC(1,J)
      N2      = NDC(2,J)
      N3      = NDC(3,J)
      N4      = NDC(4,J)
      XBAR    = THIRD*(X(1,N1)  +X(1,N2)  +X(1,N3))
      YBAR    = THIRD*(X(2,N1)  +X(2,N2)  +X(2,N3))
      ZBAR    = THIRD*(X(3,N1)  +X(3,N2)  +X(3,N3))
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      VTET1   = RNX*(X(1,N4)  -XBAR)  +RNY*(X(2,N4)  -YBAR) &
                                      +RNZ*(X(3,N4)  -ZBAR)
      XBAR    = XBAR  +THIRD*(XC(1,N1)  +XC(1,N2)  +XC(1,N3))
      YBAR    = YBAR  +THIRD*(XC(2,N1)  +XC(2,N2)  +XC(2,N3))
      ZBAR    = ZBAR  +THIRD*(XC(3,N1)  +XC(3,N2)  +XC(3,N3))
      RNX  = COFACT (X(2,N1)+XC(2,N1),X(2,N2)+XC(2,N2),X(2,N3)+XC(2,N3), &
                     X(3,N1)+XC(3,N1),X(3,N2)+XC(3,N2),X(3,N3)+XC(3,N3))
      RNY  = COFACT (X(3,N1)+XC(3,N1),X(3,N2)+XC(3,N2),X(3,N3)+XC(3,N3), &
                     X(1,N1)+XC(1,N1),X(1,N2)+XC(1,N2),X(1,N3)+XC(1,N3))
      RNZ  = COFACT (X(1,N1)+XC(1,N1),X(1,N2)+XC(1,N2),X(1,N3)+XC(1,N3), &
                     X(2,N1)+XC(2,N1),X(2,N2)+XC(2,N2),X(2,N3)+XC(2,N3))
      VTET2   = RNX*(X(1,N4)  +XC(1,N4)  -XBAR) &
               +RNY*(X(2,N4)  +XC(2,N4)  -YBAR) &
               +RNZ*(X(3,N4)  +XC(3,N4)  -ZBAR)
      STRTET  = ABS(VTET2/VTET1)
      STRMAX  = MAX(STRMAX,STRTET)
      STRMIN  = MIN(STRMIN,STRTET)
      IF (VTET1*VTET2.LE.0.) NNEG = NNEG  +1
  130 CONTINUE
!
      CALL SNGLAR (NNODE,X,XC,NCELL,NDC,NBH,SIGMN,SIGMX, &
                   SIG1,SIG2,SIG3,CNDMIN,CNDMAX)
!
      WRITE (6,700) STRMAX,STRMIN,SIGMX,SIGMN,CNDMAX,CNDMIN
      IF (NNEG.GT.0) WRITE (6,705) NNEG
      IF (NNEG.GT.0) THEN
         NFAIL = 1
         RETURN
      ENDIF
!
!
!     ADD COMPUTED DISPLACEMENTS TO THE OLD POINT POSITIONS
!
      DO 140 N=1,NNODE
      IF (ITYP(N).LT.0) GO TO 140
      IPOINT(N)   = 0
      X(1,N)      = X(1,N)  +XC(1,N)
      X(2,N)      = X(2,N)  +XC(2,N)
      X(3,N)      = X(3,N)  +XC(3,N)
  140 CONTINUE
      WRITE (6,920) L
  920 FORMAT(/'SUB-ITERATION ',I3,' IS COMPLETE'/)
  150 CONTINUE
!
!     DETERMINE MINIMUM AND MAXIMUM EXTENT OF NEW GEOMETRY
!
      NPASS      = 0
      DO 155 N=1,NNODE
      IF (ITYP(N).LT.0) GO TO 155
      IF (NPASS.EQ.0) THEN
         NPASS      = 1
         XFAR(1)    = X(1,N)
         XFAR(2)    = X(1,N)
         YFAR(1)    = X(2,N)
         YFAR(2)    = X(2,N)
         ZFAR(1)    = X(3,N)
         ZFAR(2)    = X(3,N)
      ELSE
         XFAR(1)    = MIN(XFAR(1),X(1,N))
         XFAR(2)    = MAX(XFAR(2),X(1,N))
         YFAR(1)    = MIN(YFAR(1),X(2,N))
         YFAR(2)    = MAX(YFAR(2),X(2,N))
         ZFAR(1)    = MIN(ZFAR(1),X(3,N))
         ZFAR(2)    = MAX(ZFAR(2),X(3,N))
      ENDIF
  155 CONTINUE
      WRITE (6,925) XFAR(1),XFAR(2),YFAR(1),YFAR(2),ZFAR(1),ZFAR(2)
  925 FORMAT('XFAR ',F10.2,1X,F10.2,' YFAR ',F10.2,1X,F10.2, &
             ' ZFAR ',F10.2,1X,F10.2)
!
!     COMPUTE VOLUME, CIRCUMCENTER AND CIRCUMRADIUS
!     FOR DEFORMED TETRAHEDRA
!
      DO 160 J=1,NCELL
      IF (NBH(1,J).EQ.0) GO TO 160
      N1      = NDC(1,J)
      N2      = NDC(2,J)
      N3      = NDC(3,J)
      N4      = NDC(4,J)
      IF (N4.EQ.-1) GO TO 160

      CALL CIRCUM (X,N1,N2,N3,N4,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (ISMALL.EQ.1) GO TO 310
      XCEN(J) = XCN
      YCEN(J) = YCN
      ZCEN(J) = ZCN
      VOL(J)  = VL
      RC(J)   = RADC
      AREA    = TETAR (J,X,NDC)
      RAT(J)  = RC(J)*AREA/VOL(J)
  160 CONTINUE
!      TIM         = SECOND (0)
!      WRITE (6,710) TIM
      RETURN
  310 WRITE (6,610)
      STOP
  600 FORMAT(/5X,'MAXIMUM RESIDUAL IS ',F10.6, &
                 ' MAXIMUM ITERATIONS = ',I5)
  610 FORMAT(///5X,'AT LEAST ONE NEW TETRAHEDRON HAS TOO SMALL A VOLUME' &
               /5X,'PROGRAM STOPPED IN TETMV')
  700 FORMAT(/5X,'****************************************************'/ &
              5X,'**                                                **'/ &
              5X,'**    MAXIMUM CELL STRETCHING IS  ',F10.4,'      **'/ &
              5X,'**    MAXIMUM CELL COMPRESSION IS ',F10.4,'      **'/ &
              5X,'**                                                **'/ &
              5X,'**      MAXIMUM SINGULAR VALUE IS ',F10.4,'      **'/ &
              5X,'**      MINIMUM SINGULAR VALUE IS ',F10.4,'      **'/ &
              5X,'**                                                **'/ &
              5X,'**      MAXIMUM CONDITION NUMBER  ',F10.4,'      **'/ &
              5X,'**      MINIMUM CONDITION NUMBER  ',F10.4,'      **'/ &
              5X,'**                                                **'/ &
              5X,'****************************************************'/)
  705 FORMAT(/5X,'THERE ARE ',I6,' CELLS THAT HAVE BEEN INVERTED.'/ &
              5X,'THE NEW MESH IS THEREFORE INVALID.')
!  710 FORMAT(/5X,'TIME = ',F7.1,' SECONDS'/)
      END SUBROUTINE TETMV





!
!     ******************************************************************
!
      SUBROUTINE TETMOD (X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                         NDG,IDGP,NDGP,NEDGE,NFCE,NBFACE, &
                         VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                         SIG1,SIG2,SIG3,NVCNT,RESID,COUNT,FAC, &
                         IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                         XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                         XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                         IPOINT,NPOINT,IFLAG,NFLAG, &
                         DX,DY,DZ,DS,VLT,IRING,NTETKP,NFAD,NEWC, &
                         NBHKP,IEDKP,LNBR,ISHK,MNBR,KSHK,NPP, &
                         NFILL,NEWCEL,NTRI, &
                         NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                         NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                         LISTF,VOLMIN,RCMX,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *   MODIFY THE MESH TOPOLOGY BY COARSENING AND ENRICHMENT IN     *
!     *   ORDER TO IMPROVE THE QUALITY OF THE TETRAHEDRAL ELEMENTS.    *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   2002                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,NBFACE,NCELL,NEDGE,NNODE
      INTEGER :: LNKUP(*),LNKDN(*)
      INTEGER :: NSRCH(*),KSRCH(*)
      INTEGER :: NVCNT(*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: ITYP(*),NPTET(*),LISTF(*),NACPT(*)
      INTEGER :: NPOINT(*),IPOINT(*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: NFCE(3,*)
      INTEGER :: NTRI(3,*),NFILL(*),NEWNBH(4,*),NOLD(*),NEWCEL(*), &
                 NSHAKE(*),NCAV(4,*),NEDGRM(*)
      INTEGER :: LDEL(*),NCAVFC(3,*),IKEEP(*)      
      INTEGER :: IRING(*),NTETKP(*),NFAD(3,*),NEWC(*), &
                 NBHKP(3,*),IEDKP(4,*),NPP(*), &
                 LNBR(*),ISHK(*),MNBR(*),KSHK(*)      
      DOUBLE PRECISION :: VOLMIN,TOLV,RCMX
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: SIG1(*),SIG2(*),SIG3(*)
      DOUBLE PRECISION :: RESID(*),COUNT(*),FAC(*)
      DOUBLE PRECISION :: XOCTR(2,*),YOCTR(2,*),ZOCTR(2,*), &
                          XHOLD(2,*),YHOLD(2,*),ZHOLD(2,*), &
                          XFAR(2),YFAR(2),ZFAR(2),XKEEP(2), &
                          YKEEP(2),ZKEEP(2)
      DOUBLE PRECISION :: XC(*),YC(*),ZC(*),V(*),RAD(*),RCRIN(*)
      DOUBLE PRECISION :: DX(*),DY(*),DZ(*),DS(*),VLT(*)
      
      INTEGER :: JFIRST,JLAST,LC,NBPTS,NTRACK
!
!     ******************************************************************
!
!     COMPUTE DENSITY FUNCTION AT EACH MESH POINT.
!
      CALL DENSFN (X,NNODE,NFCE,NBFACE,NEDGE,NDG,DENS,RESID,COUNT)
!
!     COMPUTE EDGE VALENCE OF EACH MESH POINT.
! 
      CALL EDGLEN (X,NNODE,ITYP,NEDGE,NDG,DENS,NVCNT)
!
!     COARSEN MESH BY EDGE COLLAPSE.
!
      DO 25 LC=1,3
!
      CALL COARSN (LC,X,NNODE,NDC,NBH,IPROT,NCELL, &
                   ITYP,XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                   NVCNT,DENS,IFLAG,NFLAG,NPTET, &
                   NEDGE,NDG,IDGP,NDGP, &
                   NOCTR,IOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE,NREF, &
                   FAC,SIG1,SIG2,SIG3, &
                   KSRCH,NSRCH,IRING,NTETKP, &
                   LNBR,ISHK,MNBR,KSHK,TOLV)
!
!     OPTIMIZE MESH QUALITY BY FACE/EDGE SWAPPING.
!
      CALL SMOOTH (X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                   NDG,IDGP,NDGP,NEDGE, &
                   VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                   IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                   XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                   XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                   IPOINT,NPOINT,IFLAG,NFLAG, &
                   DX,DY,DZ,DS,VLT,IRING,NTETKP,NFAD,NEWC, &
                   NBHKP,IEDKP,LNBR,ISHK,MNBR,KSHK,NPP, &
                   NFILL,NEWCEL,NTRI, &
                   NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                   NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                   LISTF,VOLMIN,RCMX,TOLV)
!
!     RE-COMPUTE EDGE VALENCE OF EACH MESH POINT.
! 
      CALL EDGLEN (X,NNODE,ITYP,NEDGE,NDG,DENS,NVCNT)

   25 CONTINUE
!
!     RE-COMPUTE DENSITY FUNCTION AT EACH MESH POINT.
!
      CALL DENSFN(X,NNODE,NFCE,NBFACE,NEDGE,NDG,DENS,RESID,COUNT)
!
!     ENRICH MESH INTERIOR BASED ON A COMPARISON BETWEEN THE
!     ACTUAL POINT DENSITY AND THE VALUE OF THE POINT 
!     DENSITY FUNCTION.
!
      CALL VOLPUT (X,ITYP,NBPTS,NNODE,NDC,NBH,IPROT,NCELL, &
                   NDG,IDGP,NDGP,NEDGE, &
                   VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                   IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                   XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                   XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                   IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                   NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                   NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                   JLAST,JFIRST,NTRACK,VOLMIN,RCMX,TOLV) 
!
!     OPTIMIZE MESH QUALITY BY FACE/EDGE SWAPPING.
!
      CALL SMOOTH (X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                   NDG,IDGP,NDGP,NEDGE, &
                   VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                   IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                   XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                   XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                   IPOINT,NPOINT,IFLAG,NFLAG, &
                   DX,DY,DZ,DS,VLT,IRING,NTETKP,NFAD,NEWC, &
                   NBHKP,IEDKP,LNBR,ISHK,MNBR,KSHK,NPP, &
                   NFILL,NEWCEL,NTRI, &
                   NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                   NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                   LISTF,VOLMIN,RCMX,TOLV)
!
!     RE-COMPUTE EDGE VALENCE OF EACH MESH POINT.
! 
      CALL EDGLEN (X,NNODE,ITYP,NEDGE,NDG,DENS,NVCNT)

      RETURN
      END SUBROUTINE TETMOD





!
!     ******************************************************************
!
      SUBROUTINE OCTFIL (N,X,NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR)
!
!     ******************************************************************
!     *                                                                *
!     *   INSERT POINT N INTO OCTREE STRUCTURE                         *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *   COPYRIGHT (C) TIM BAKER     1994                             *
!     *                                                                *
!     ******************************************************************
!
!
      IMPLICIT NONE

      INTEGER :: IOCTR,N
      INTEGER :: NREF(*),NLINK(*),NOCTR(2,*)
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)
      DOUBLE PRECISION :: XROOT(2),YROOT(2),ZROOT(2)
      
      INTEGER :: I,IROOT,J,K,L,LTEST,NEXT,NTEST,NXSGN,NYSGN,NZSGN
      INTEGER :: JJ(8),NSTORE(8)      
      DOUBLE PRECISION :: TOL,TOLPT, &
                          XHALF,XHIGH,XLOW,XSGN,XSHIFT,XSIZE, & 
                          YHALF,YHIGH,YLOW,YSGN,YSHIFT,YSIZE, &
                          ZHALF,ZHIGH,ZLOW,ZSGN,ZSHIFT,ZSIZE
!
!     ******************************************************************
!
      TOLPT      = 1.E-6
      TOL        = 1.000001
      IF (X(1,N).LT.XFAR(1)-TOLPT.OR.X(1,N).GT.XFAR(2)+TOLPT) GO TO 200
      IF (X(2,N).LT.YFAR(1)-TOLPT.OR.X(2,N).GT.YFAR(2)+TOLPT) GO TO 200
      IF (X(3,N).LT.ZFAR(1)-TOLPT.OR.X(3,N).GT.ZFAR(2)+TOLPT) GO TO 200
   5  I          = 1
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
      NREF(N)    = I
      IF (L.GT.1) GO TO 15
      NOCTR(1,I) = N
      RETURN
   15 NEXT       = NOCTR(1,I)
   20 NTEST      = NLINK(NEXT)
      IF (NTEST.EQ.0) GO TO 25
      NEXT       = NTEST
      GO TO 20
   25 NLINK(NEXT) = N
      RETURN
!
!     OCTANT NOW CONTAINS NINE POINTS. FORM EIGHT NEW OCTANTS
!     AND AMEND TREE STRUCTURE.
!
   30 DO 35 K=1,8
      IOCTR      = IOCTR  +1
      IF (IOCTR.GT.MXOCTR) GO TO 210
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
      XSHIFT     = X(1,J)  -.5*(XLOW  +XHIGH)
      XSIZE      = MAX(1.0D0-9,ABS(XSHIFT))
      NXSGN      = (INT(REAL(TOL*XSHIFT/XSIZE))  +1)/2*2
!     NXSGN      = (IFIX(REAL(TOL*XSHIFT/XSIZE))  +1)/2*2
      YSHIFT     = X(2,J)  -.5*(YLOW  +YHIGH)
      YSIZE      = MAX(1.0D0-9,ABS(YSHIFT))
      NYSGN      = (INT(REAL(TOL*YSHIFT/YSIZE))  +1)/2*2
!     NYSGN      = (IFIX(REAL(TOL*YSHIFT/YSIZE))  +1)/2*2
      ZSHIFT     = X(3,J)  -.5*(ZLOW  +ZHIGH)
      ZSIZE      = MAX(1.0D0-9,ABS(ZSHIFT))
      NZSGN      = (INT(REAL(TOL*ZSHIFT/ZSIZE))  +1)/2*2
!     NZSGN      = (IFIX(REAL(TOL*ZSHIFT/ZSIZE))  +1)/2*2
      L          = 1  +NXSGN/2  +NYSGN  +2*NZSGN
      JJ(L)      = JJ(L)  +1
      NLINK(J)   = 0
      NREF(J)    = IOCTR  -8  +L
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
      XHALF      = .5*(XLOW  +XHIGH)
      XSGN       = MOD(LTEST,2)
      XLOW       = XSGN*XLOW   +(1.  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1.  -XSGN)*XHIGH
      YHALF      = .5*(YLOW  +YHIGH)
      YSGN       = ISIGN(1,2*MOD(LTEST-1,4)-3)
      YLOW       = .5*((1.  -YSGN)*YLOW   +(1.  +YSGN)*YHALF)
      YHIGH      = .5*((1.  -YSGN)*YHALF  +(1.  +YSGN)*YHIGH)
      ZHALF      = .5*(ZLOW  +ZHIGH)
      ZSGN       = ISIGN(1,2*LTEST-9)
      ZLOW       = .5*((1.  -ZSGN)*ZLOW   +(1.  +ZSGN)*ZHALF)
      ZHIGH      = .5*((1.  -ZSGN)*ZHALF  +(1.  +ZSGN)*ZHIGH)
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
   80 XHALF      = .5*(XLOW  +XHIGH)
      YHALF      = .5*(YLOW  +YHIGH)
      ZHALF      = .5*(ZLOW  +ZHIGH)
      XSHIFT     = X(1,N)  -XHALF
      XSIZE      = MAX(1.0D0-9,ABS(XSHIFT))
      NXSGN      = (INT(REAL(TOL*XSHIFT/XSIZE))  +1)/2*2
!     NXSGN      = (IFIX(REAL(TOL*XSHIFT/XSIZE))  +1)/2*2
      YSHIFT     = X(2,N)  -YHALF
      YSIZE      = MAX(1.0D0-9,ABS(YSHIFT))
      NYSGN      = (INT(REAL(TOL*YSHIFT/YSIZE))  +1)/2*2
!     NYSGN      = (IFIX(REAL(TOL*YSHIFT/YSIZE))  +1)/2*2
      ZSHIFT     = X(3,N)  -ZHALF
      ZSIZE      = MAX(1.0D0-9,ABS(ZSHIFT))
      NZSGN      = (INT(REAL(TOL*ZSHIFT/ZSIZE))  +1)/2*2
!     NZSGN      = (IFIX(REAL(TOL*ZSHIFT/ZSIZE))  +1)/2*2
      L          = 1  +NXSGN/2  +NYSGN  +2*NZSGN
      I          = -NOCTR(1,I)  +L  -1
      XSGN       = MOD(L,2)
      XLOW       = XSGN*XLOW   +(1.  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1.  -XSGN)*XHIGH
      YSGN       = ISIGN(1,2*MOD(L-1,4)-3)
      YLOW       = .5*((1.  -YSGN)*YLOW   +(1.  +YSGN)*YHALF)
      YHIGH      = .5*((1.  -YSGN)*YHALF  +(1.  +YSGN)*YHIGH)
      ZSGN       = ISIGN(1,2*L-9)
      ZLOW       = .5*((1.  -ZSGN)*ZLOW   +(1.  +ZSGN)*ZHALF)
      ZHIGH      = .5*((1.  -ZSGN)*ZHALF  +(1.  +ZSGN)*ZHIGH)
      IF (NOCTR(1,I).LT.0) GO TO 80
      L          = 0
      IF (NOCTR(1,I).EQ.0) GO TO 10
      NEXT       = NOCTR(1,I)
   85 NEXT       = NLINK(NEXT)
      L          = L  +1
      IF (NEXT.NE.0) GO TO 85
      GO TO 10
  200 WRITE (6,600) N

      WRITE (6,900) N,X(1,N),X(2,N),X(3,N), &
                    XFAR(1),XFAR(2),YFAR(1),YFAR(2),ZFAR(1),ZFAR(2)
  900 FORMAT('N ',I4,' X ',F8.3,' Y ',F8.3,' Z ',F8.3/ &
             'XFAR',2(1X,F8.3),' YFAR',2(1X,F8.3),' ZFAR',2(1X,F8.3))

      STOP
  210 WRITE (6,610)
      STOP
  600 FORMAT(5X,'POINT N= ',I6,'      LIES OUTSIDE CONVEX HULL')
  610 FORMAT(5X,'DIMENSION OF NOCTR ARRAY EXCEEDED'/ &
             5X,'INCREASE SIZE OF MOCTR')
      END SUBROUTINE OCTFIL





!
!     ******************************************************************
!
      SUBROUTINE OCTFND (NCODE,NCLOSE,XPT,YPT,ZPT,N1,N2,N3,DISMIN, &
                         X,NOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE, &
                         XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                         XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH)
!
!     ******************************************************************
!     *                                                                *
!     *  FIND NEAREST POINT NCLOSE IN OCTREE TO (XPT,YPT,ZPT)          *
!     *  THAT LIES IN FRONT OF TRIANGLE (N1,N2,N3)
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *   COPYRIGHT (C) TIM BAKER     1994                             *
!     *                                                                *
!     ******************************************************************
!
!
      IMPLICIT NONE

      INTEGER :: NCODE,NCLOSE,N1,N2,N3
      INTEGER :: IDONE(*),NLINK(*),NOCTR(2,*)
      INTEGER :: NSRCH(*),KSRCH(*)
      DOUBLE PRECISION :: DISMIN,XPT,YPT,ZPT
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: XOCTR(2,*),YOCTR(2,*),ZOCTR(2,*), &
                          XHOLD(2,*),YHOLD(2,*),ZHOLD(2,*), &
                          XKEEP(2),YKEEP(2),ZKEEP(2),XFAR(2), &
                          YFAR(2),ZFAR(2)

      INTEGER :: I,IC,IFLAG,II,ITRY,J,JJ,K,KC,L,LFLAG,NEXT,NXSGN, &
                 NYSGN,NZSGN
      DOUBLE PRECISION :: C1,C2,C3,DET,DIST,DMIN,TOL,TOLPT, &
                          XHALF,XHIGH,XL,XLOW,XSHIFT,XSIZE,XSGN,XU, &
                          YHALF,YHIGH,YL,YLOW,YSHIFT,YSIZE,YSGN,YU, &
                          ZHALF,ZHIGH,ZL,ZLOW,ZSHIFT,ZSIZE,ZSGN,ZU
!
!     ******************************************************************
!
      TOLPT      = 1.E-6
      TOL        = 1.000001
      IF (XPT.LT.XFAR(1)-TOLPT.OR.XPT.GT.XFAR(2)+TOLPT) GO TO 200
      IF (YPT.LT.YFAR(1)-TOLPT.OR.YPT.GT.YFAR(2)+TOLPT) GO TO 200
      IF (ZPT.LT.ZFAR(1)-TOLPT.OR.ZPT.GT.ZFAR(2)+TOLPT) GO TO 200
      IF (N1.EQ.0) GO TO 5
!
!     FORM INWARD POINTING NORMAL FOR TRIANGLE (N1,N2,N3)
!
      C1      = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      C2      = COFACT (X(1,N1),X(1,N2),X(1,N3),X(3,N1),X(3,N2),X(3,N3))
      C3      = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
!
!     FIRST FIND OCTANT WHICH CONTAINS POINT (XPT,YPT,ZPT)
!
   5  I          = 1
      XLOW       = XFAR(1)
      XHIGH      = XFAR(2)
      YLOW       = YFAR(1)
      YHIGH      = YFAR(2)
      ZLOW       = ZFAR(1)
      ZHIGH      = ZFAR(2)
      IF (NOCTR(1,I).LT.0) GO TO 20
      NEXT       = NOCTR(1,I)
   10 NEXT       = NLINK(NEXT)
      IF (NEXT.NE.0) GO TO 10
      GO TO 50
!
!     LOCATE SUB-OCTANT IN WHICH POINT LIES
!
   20 XHALF      = .5*(XLOW  +XHIGH)
      YHALF      = .5*(YLOW  +YHIGH)
      ZHALF      = .5*(ZLOW  +ZHIGH)
      XSHIFT     = XPT  -XHALF
      XSIZE      = MAX(1.0D0-9,ABS(XSHIFT))
      NXSGN      = (INT(REAL(TOL*XSHIFT/XSIZE))  +1)/2*2
!     NXSGN      = (IFIX(REAL(TOL*XSHIFT/XSIZE))  +1)/2*2
      YSHIFT     = YPT  -YHALF
      YSIZE      = MAX(1.0D0-9,ABS(YSHIFT))
      NYSGN      = (INT(REAL(TOL*YSHIFT/YSIZE))  +1)/2*2
!     NYSGN      = (IFIX(REAL(TOL*YSHIFT/YSIZE))  +1)/2*2
      ZSHIFT     = ZPT  -ZHALF
      ZSIZE      = MAX(1.0D0-9,ABS(ZSHIFT))
      NZSGN      = (INT(REAL(TOL*ZSHIFT/ZSIZE))  +1)/2*2
!     NZSGN      = (IFIX(REAL(TOL*ZSHIFT/ZSIZE))  +1)/2*2
      L          = 1  +NXSGN/2  +NYSGN  +2*NZSGN
      I          = -NOCTR(1,I)  +L  -1
      XSGN       = MOD(L,2)
      XLOW       = XSGN*XLOW   +(1.  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1.  -XSGN)*XHIGH
      YSGN       = ISIGN(1,2*MOD(L-1,4)-3)
      YLOW       = .5*((1.  -YSGN)*YLOW   +(1.  +YSGN)*YHALF)
      YHIGH      = .5*((1.  -YSGN)*YHALF  +(1.  +YSGN)*YHIGH)
      ZSGN       = ISIGN(1,2*L-9)
      ZLOW       = .5*((1.  -ZSGN)*ZLOW   +(1.  +ZSGN)*ZHALF)
      ZHIGH      = .5*((1.  -ZSGN)*ZHALF  +(1.  +ZSGN)*ZHIGH)
      IF (NOCTR(1,I).LT.0) GO TO 20
      IF (NOCTR(1,I).EQ.0) GO TO 50
      NEXT       = NOCTR(1,I)
   25 NEXT       = NLINK(NEXT)
      IF (NEXT.NE.0) GO TO 25
!
!     SEARCH FOR POINT NCLOSE THAT IS NEAREST TO (XPT,YPT,ZPT)
!
   50 NCLOSE      = 1
      DISMIN      = (XFAR(2)  -XFAR(1))**2  +(YFAR(2)  -YFAR(1))**2 &
                   +(ZFAR(2)  -ZFAR(1))**2
      IF (NOCTR(1,I).EQ.0) GO TO 65
      J           = NOCTR(1,I)
   58 IF (IDONE(J).EQ.0.AND.NCODE.EQ.0) GO TO 60
      IF (IDONE(J).EQ.2) GO TO 60
      IF (N1.EQ.0) GO TO 59
      DET       = C1*(X(1,J)  -XPT)  -C2*(X(2,J)  -YPT) &
                                     +C3*(X(3,J)  -ZPT)
      IF (DET.LT.TOLPT) GO TO 60
   59 DIST        = (XPT  -X(1,J))**2  +(YPT  -X(2,J))**2 &
                   +(ZPT  -X(3,J))**2
      IF (DIST.GE.DISMIN) GO TO 60
      NCLOSE      = J
      DISMIN      = DIST
   60 J           = NLINK(J)
      IF (J.NE.0) GO TO 58
   65 DMIN        = SQRT(DISMIN)
      XL          = XPT  -DMIN
      XU          = XPT  +DMIN
      YL          = YPT  -DMIN
      YU          = YPT  +DMIN
      ZL          = ZPT  -DMIN
      ZU          = ZPT  +DMIN
   70 IF (XL.GT.XLOW.AND.XU.LT.XHIGH.AND. &
          YL.GT.YLOW.AND.YU.LT.YHIGH.AND. &
          ZL.GT.ZLOW.AND.ZU.LT.ZHIGH) RETURN
      IF (I.EQ.1) RETURN
      IFLAG       = I
      I           = NOCTR(2,I)
      LFLAG       = IFLAG  +NOCTR(1,I)  +1
      XSGN        = MOD(LFLAG,2)
      XOCTR(1,1)  = (2.  -XSGN)*XLOW  -(1.  -XSGN)*XHIGH
      XOCTR(2,1)  = (1.  +XSGN)*XHIGH  -XSGN*XLOW
      YSGN        = ISIGN(1,2*MOD(LFLAG-1,4)-3)
      YOCTR(1,1)  = .5*((3.  +YSGN)*YLOW  -(1.  +YSGN)*YHIGH)
      YOCTR(2,1)  = .5*((3.  -YSGN)*YHIGH  -(1.  -YSGN)*YLOW)
      ZSGN        = ISIGN(1,2*LFLAG-9)
      ZOCTR(1,1)  = .5*((3.  +ZSGN)*ZLOW  -(1.  +ZSGN)*ZHIGH)
      ZOCTR(2,1)  = .5*((3.  -ZSGN)*ZHIGH  -(1.  -ZSGN)*ZLOW)
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
      XHALF      = .5*(XOCTR(1,J)  +XOCTR(2,J))
      XSGN       = MOD(K,2)
      XLOW       = XSGN*XOCTR(1,J)   +(1.  -XSGN)*XHALF
      XHIGH      = XSGN*XHALF  +(1.  -XSGN)*XOCTR(2,J)
      YHALF      = .5*(YOCTR(1,J)  +YOCTR(2,J))
      YSGN       = ISIGN(1,2*MOD(K-1,4)-3)
      YLOW       = .5*((1.  -YSGN)*YOCTR(1,J)   +(1.  +YSGN)*YHALF)
      YHIGH      = .5*((1.  -YSGN)*YHALF  +(1.  +YSGN)*YOCTR(2,J))
      ZHALF      = .5*(ZOCTR(1,J)  +ZOCTR(2,J))
      ZSGN       = ISIGN(1,2*K-9)
      ZLOW       = .5*((1.  -ZSGN)*ZOCTR(1,J)   +(1.  +ZSGN)*ZHALF)
      ZHIGH      = .5*((1.  -ZSGN)*ZHALF  +(1.  +ZSGN)*ZOCTR(2,J))
      IF (XL.GT.XHIGH.OR.XU.LT.XLOW.OR. &
          YL.GT.YHIGH.OR.YU.LT.YLOW.OR. &
          ZL.GT.ZHIGH.OR.ZU.LT.ZLOW) GO TO 90
      IF (NOCTR(1,ITRY).GE.0) GO TO 80
      IC          = IC  +1
      IF (IC.GT.MXTEST) GO TO 210
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
   82 IF (IDONE(JJ).EQ.0.AND.NCODE.EQ.0) GO TO 85
      IF (IDONE(JJ).EQ.2) GO TO 85
      IF (N1.EQ.0) GO TO 83
      DET       = C1*(X(1,JJ)  -XPT)  -C2*(X(2,JJ)  -YPT) &
                                      +C3*(X(3,JJ)  -ZPT)
      IF (DET.LT.TOLPT) GO TO 85
   83 DIST        = (XPT  -X(1,JJ))**2  +(YPT  -X(2,JJ))**2 &
                   +(ZPT  -X(3,JJ))**2
      IF (DIST.GE.DISMIN) GO TO 85
      NCLOSE      = JJ
      DISMIN      = DIST
   85 JJ          = NLINK(JJ)
      IF (JJ.NE.0) GO TO 82
      DMIN        = SQRT(DISMIN)
      XL          = XPT  -DMIN
      XU          = XPT  +DMIN
      YL          = YPT  -DMIN
      YU          = YPT  +DMIN
      ZL          = ZPT  -DMIN
      ZU          = ZPT  +DMIN
   90 CONTINUE
      IF (IC.GT.0) GO TO 92
      XLOW        = XKEEP(1)
      XHIGH       = XKEEP(2)
      YLOW        = YKEEP(1)
      YHIGH       = YKEEP(2)
      ZLOW        = ZKEEP(1)
      ZHIGH       = ZKEEP(2)
      GO TO 70
   92 KC          = IC
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
             5X,'LIES OUTSIDE CONVEX HULL')
  610 FORMAT(5X,'DIMENSION OF NSRCH AND KSRCH ARRAYS EXCEEDED'/ &
             5X,'IN ROUTINE OCTFND. INCREASE SIZE OF MTEST')
      END SUBROUTINE OCTFND




!
!     ******************************************************************
!
      SUBROUTINE OCTRMV (N,X,NOCTR,NLINK,IDONE,NREF)
!
!     ******************************************************************
!     *                                                                *
!     *  REMOVE POINT N FROM OCTREE DATA STRUCTURE                     *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *   COPYRIGHT (C) TIM BAKER     1998                             *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: N
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      DOUBLE PRECISION :: X(3,*)
      
      INTEGER :: I,NEXT,NPRE
!
!     ******************************************************************
!
      I          = NREF(N)
      NREF(N)    = 0
      IF (I.EQ.0) GO TO 200
      NEXT       = NOCTR(1,I)
      IF (NEXT.NE.N) GO TO 10
      NEXT       = NLINK(N)
      NOCTR(1,I) = NEXT
      NLINK(N)   = 0
      IDONE(N)   = 0
      RETURN
   10 IF (NEXT.EQ.0) GO TO 210
      NPRE       = NEXT
      NEXT       = NLINK(NEXT)
      IF (NEXT.NE.N) GO TO 10
      NLINK(NPRE) = NLINK(N)
      NLINK(N)  = 0
      IDONE(N)   = 0
      RETURN
  200 WRITE (6,600) N,X(1,N),X(2,N),X(3,N)
      STOP
  210 WRITE (6,610) N,X(1,N),X(2,N),X(3,N)
      STOP
  600 FORMAT(5X,'POINT ',I6,' HAS AN NREF VALUE OF ZERO'/ &
             5X,'X = ',F10.4,' Y = ',F10.4,' Z = ',F10.4)
  610 FORMAT(5X,'UNABLE TO FIND POINT WITH ADDRESS ',I6,' IN OCTREE'/ &
             5X,'X = ',F10.4,' Y = ',F10.4,' Z = ',F10.4)
      END SUBROUTINE OCTRMV




!
!     ******************************************************************
!
      SUBROUTINE EDGERM (NEDG,NDG,IDGP,NDGP)
!
!     ******************************************************************
!     *                                                                *
!     *   REMOVE EDGE NEDG FROM LINKED LIST OF EDGES                   *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      INTEGER :: NEDG
      INTEGER :: NDG(2,*),IDGP(*),NDGP(*)
      
      INTEGER :: I,IEDG,N1
!
!     ******************************************************************
!
      IF (NEDG.EQ.0) GO TO 200
      N1         = MIN(NDG(1,NEDG),NDG(2,NEDG))
      NDG(1,NEDG) = -1
      NDG(2,NEDG) = -1
      I          = IDGP(N1)
      IF (I.NE.NEDG) GO TO 10
      I          = NDGP(I)
      IDGP(N1)   = I
      IF (I.EQ.0) RETURN
      IEDG       = I
      GO TO 20
   10 IEDG       = I
      I          = NDGP(I)
      IF (I.EQ.0) GO TO 210
      IF (I.NE.NEDG) GO TO 10
   20 I          = NDGP(I)
      NDGP(IEDG) = I
      IF (I.EQ.0) RETURN
      IEDG       = I
      GO TO 20
  200 WRITE (6,600)
      STOP
  210 WRITE (6,610)
      STOP
  600 FORMAT(//5X,'EDGE REQUESTED HAS ADDRESS ZERO'/ &
               5X,'PROGRAM STOPPED IN EDGERM')
  610 FORMAT(//5X,'UNABLE TO FIND EDGE IN LINKED LIST'/ &
               5X,'PROGRAM STOPPED IN EDGERM')
      END SUBROUTINE EDGERM




!
!     ******************************************************************
!
      SUBROUTINE TREE (PROP,NLEFT,NRIGHT,NBACK,LISTF,NACPT,NTOT, &
                       NBH,IPROT,NCELL)
!
!     ******************************************************************
!     *                                                                *
!     *  CREATE A BINARY TREE AND ASSEMBLE AN ORDERED LIST WITH        *
!     *  RESPECT TO THE SIZE OF THE SCALAR PROPERTY PROP               *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *  COPYRIGHT (C) TIM BAKER    1994                               *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      INTEGER :: NCELL,NTOT
      INTEGER :: NBH(4,*),IPROT(*)
      INTEGER :: NLEFT(*),NRIGHT(*),NBACK(*),LISTF(*),NACPT(*)
      DOUBLE PRECISION :: PROP(*)
      
      INTEGER :: J,JBACK,JJ,JNEXT,JPROBE,JSTART,L,NCNT,NJ
      DOUBLE PRECISION :: TOL
!
!     ******************************************************************
!
      TOL       = 1.E-12
!
!     CREATE BINARY TREE
!
      NTOT      = 0
      J         = 0
   10 J         = J  +1
      IF (J.GT.NCELL) GO TO 90
      IF (IPROT(J).EQ.1) GO TO 10
      IF (NBH(1,J).EQ.0) GO TO 10
      JSTART    = J
      NBACK(J)  = 0
      NLEFT(J)  = 0
      NRIGHT(J) = 0
      NTOT      = 1
   20 J         = J  +1
      IF (J.GT.NCELL) GO TO 90
      IF (IPROT(J).EQ.1) GO TO 20
      IF (NBH(1,J).EQ.0) GO TO 20
      JJ        = JSTART
   30 JNEXT     = NLEFT(JJ)
      IF (PROP(J).GT.PROP(JJ)) JNEXT = NRIGHT(JJ)
      IF (JNEXT.EQ.0) GO TO 40
      JJ        = JNEXT
      GO TO 30
   40 IF (PROP(J).GT.PROP(JJ)) GO TO 50
      NLEFT(JJ) = J
      NBACK(J)  = -JJ
      GO TO 60
   50 NRIGHT(JJ) = J
      NBACK(J)   = JJ
   60 NLEFT(J)   = 0
      NRIGHT(J)  = 0
      NTOT      = NTOT  +1
      GO TO 20
!
!     ASSEMBLE ORDERED LIST
!
   90 IF (NTOT.EQ.0) RETURN
      NCNT      = 0
      JPROBE    = JSTART
  100 J         = NRIGHT(JPROBE)
      IF (J.EQ.0) GO TO 110
      JPROBE    = J
      GO TO 100
  110 J         = JPROBE
      NCNT      = NCNT  +1
      LISTF(NCNT) = J
      NACPT(J)  = 0
      GO TO 120
  115 NCNT      = NCNT  +1
      LISTF(NCNT) = J
      NACPT(J)  = 0
  120 JPROBE    = NLEFT(J)
      IF (JPROBE.GT.0) GO TO 100
  125 JBACK     = NBACK(J)
      IF (JBACK.EQ.0) GO TO 140
      NJ        = ISIGN(1,JBACK)
      J         = IABS(JBACK)
      IF (NJ.GT.0) GO TO 115
      GO TO 125
  140 DO 150 L=1,NCELL
      NBACK(L)  = 0
  150 CONTINUE
      RETURN
      END SUBROUTINE TREE





!
!     ******************************************************************
!
      SUBROUTINE SNGLAR (NNODE,X,XD,NCELL,NDC,NBH,SIGMN,SIGMX, &
                         SIG1,SIG2,SIG3,CNDMIN,CNDMAX)
!
!     ******************************************************************
!     *                                                                *
!     *  DETERMINE THE SINGULAR VALUES SIG1(L), SIG2(L), SIG3(L)  OF   *
!     *  THE MATRIX ASSOCIATED WITH THE DEFORMATION OF EACH            *
!     *  TETRAHEDRON L IN THE DISPLACED MESH.                          *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    X(3,N) contains the original x, y, z coordinates of the     *
!     *           Nth point in the mesh.                               *
!     *                                                                *
!     *    XD(3,N) contains the x, y, z displacements of the Nth point *
!     *            in the mesh.                                        *
!     *                                                                *
!     *    NDC(4,L) contains the vertex addresses of the 4 vertices    *
!     *             of tetrahedron L.                                  *
!     *                                                                *
!     *    NBH(4,L) contains the addresses of the 4 tetrahedra that    *
!     *             are neighbors to tetrahedron L.                    *
!     *             (NBH is required if SNGLAR is called from the      *
!     *             mesh modification routines. If SNGLAR is called    *
!     *             from CRUNCH then NBH is not needed and the IF      *
!     *             statement, in which NBH appears, can be commented  *
!     *             out.)                                              *
!     *                                                                *
!     *    NNODE    is the total number of mesh points.                *
!     *                                                                *
!     *    NCELL    is the total number of mesh tetrahedra.            *
!     *                                                                *
!     *    SIGMN and SIGMX are scalar values returned by SNGLAR which  *
!     *    which give the minimum and maximum, over the whole mesh, of *
!     *    the singular values of the deformation matrices.            *
!     *                                                                *
!     *    SIG1,SIG2,SIG3 are vector quantities returned by SNGLAR     *
!     *    which contain the singular values of the deformation matrix *
!     *    for each tetrahedron in the mesh. SIG1,SIG2 and SIG3 should *
!     *    each be given dimension MCELL.                              *
!     *    If, for any given tetrahedron,  SIG1 = 1, SIG2 = 1 and      *
!     *    SIG3 = 1 then the tetrahedron has not changed shape,        *
!     *    although it may have undergone a rotation or an isotropic   *
!     *    expansion or compression.                                   * 
!     *    If any of the singular values is close to zero then the     *
!     *    tetrahedron has become highly compressed.                   *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   2001                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      INTEGER :: NCELL,NNODE
      INTEGER  NDC(4,*),NBH(4,*)
      DOUBLE PRECISION :: CNDMAX,CNDMIN,SIGMN,SIGMX
      DOUBLE PRECISION :: X(3,*),XD(3,*),SIG1(*),SIG2(*),SIG3(*)
      
      INTEGER :: L,LMAX,LMIN,NPASS,N1,N2,N3,N4
      DOUBLE PRECISION :: A0,A1,A11,A12,A13,A2,A21,A22,A23,A31,A32,A33, &
                          B1,B2,B3,B4,B4SQ,B5,B5SQ,B6,B6SQ,COND,CT, &
                          D,DELX1,DELX2,DELX3,DELY1,DELY2,DELY3, &
                          DELZ1,DELZ2,DELZ3,DET,DSQ,DX1,DX2,DX3, &
                          DY1,DY2,DY3,DZ1,DZ2,DZ3,FAC,P1,P2,Q,QA, &
                          QROOT,R,S,SIGLOW,SIGLRG,SLIDE,SLIDCB,SLIDSQ, &
                          ST,S1X,S1Y,S1Z,S2X,S2Y,S2Z,S3X,S3Y,S3Z, &
                          THETA,THIRD,TOL
!
!     ******************************************************************
!
      TOL         = 1.E-12
      THIRD       = 1./3.
      NPASS       = 0
      DO 100 L=1,NCELL
      IF (NBH(1,L).EQ.0) GO TO 100
      IF (NDC(4,L).LE.0) GO TO 100
      N1          = NDC(1,L)        
      N2          = NDC(2,L)        
      N3          = NDC(3,L)        
      N4          = NDC(4,L)        
!
!     COMPUTE ENTRIES, COFACTORS AND DETERMINANT 
!     OF THE ORIGINAL EDGE MATRIX T
!
      DX1         = X(1,N2)  -X(1,N1)
      DY1         = X(2,N2)  -X(2,N1)
      DZ1         = X(3,N2)  -X(3,N1)
      DX2         = X(1,N3)  -X(1,N1)
      DY2         = X(2,N3)  -X(2,N1)
      DZ2         = X(3,N3)  -X(3,N1)
      DX3         = X(1,N4)  -X(1,N1)
      DY3         = X(2,N4)  -X(2,N1)
      DZ3         = X(3,N4)  -X(3,N1)
      S1X         = DY2*DZ3  -DY3*DZ2
      S1Y         = DZ2*DX3  -DZ3*DX2
      S1Z         = DX2*DY3  -DX3*DY2
      S2X         = DY3*DZ1  -DY1*DZ3
      S2Y         = DZ3*DX1  -DZ1*DX3
      S2Z         = DX3*DY1  -DX1*DY3
      S3X         = DY1*DZ2  -DY2*DZ1
      S3Y         = DZ1*DX2  -DZ2*DX1
      S3Z         = DX1*DY2  -DX2*DY1
      DET         = DX1*S1X  +DY1*S1Y  +DZ1*S1Z
      IF (ABS(DET).LT.TOL) GO TO 300
      FAC         = 1./DET
!
!     COMPUTE EDGE DISPLACEMENTS AND ENTRIES
!     FOR THE DEFORMATION MATRIX A = Tnew * T^(-1)
!
      DELX1       = XD(1,N2)  -XD(1,N1)
      DELY1       = XD(2,N2)  -XD(2,N1)
      DELZ1       = XD(3,N2)  -XD(3,N1)
      DELX2       = XD(1,N3)  -XD(1,N1)
      DELY2       = XD(2,N3)  -XD(2,N1)
      DELZ2       = XD(3,N3)  -XD(3,N1)
      DELX3       = XD(1,N4)  -XD(1,N1)
      DELY3       = XD(2,N4)  -XD(2,N1)
      DELZ3       = XD(3,N4)  -XD(3,N1)
      A11         = (DELX1*S1X  +DELX2*S2X  +DELX3*S3X)*FAC
      A12         = (DELX1*S1Y  +DELX2*S2Y  +DELX3*S3Y)*FAC
      A13         = (DELX1*S1Z  +DELX2*S2Z  +DELX3*S3Z)*FAC
      A21         = (DELY1*S1X  +DELY2*S2X  +DELY3*S3X)*FAC
      A22         = (DELY1*S1Y  +DELY2*S2Y  +DELY3*S3Y)*FAC
      A23         = (DELY1*S1Z  +DELY2*S2Z  +DELY3*S3Z)*FAC
      A31         = (DELZ1*S1X  +DELZ2*S2X  +DELZ3*S3X)*FAC
      A32         = (DELZ1*S1Y  +DELZ2*S2Y  +DELZ3*S3Y)*FAC
      A33         = (DELZ1*S1Z  +DELZ2*S2Z  +DELZ3*S3Z)*FAC
!
!     COMPUTE THE ENTRIES FOR THE POSITIVE DEFINITE MATRIX A^T * A
!
      B1          = (1.  +A11)**2  +A21*A21  +A31*A31
      B2          = A12*A12  +(1.  +A22)**2  +A32*A32
      B3          = A13*A13  +A23*A23  +(1.  +A33)**2
      B4          = A12*(1.  +A11)  +A21*(1.  +A22)  +A31*A32
      B5          = A13*(1.  +A11)  +A21*A23  +A31*(1.  +A33)
      B6          = A12*A13  +A23*(1.  +A22)  +A32*(1.  +A33)
!
!     COMPUTE THE COEFFICIENTS OF CHARACTERISTIC POLYNOMIAL
!       LAMBDA^3 + A2*LAMBDA^2 + A1*LAMBDA + A0 = 0
!
      B4SQ        = B4*B4
      B5SQ        = B5*B5
      B6SQ        = B6*B6
      A0          = B1*B6SQ  +B2*B5SQ  +B3*B4SQ  -B1*B2*B3  -2.*B4*B5*B6
      A1          = B1*B2  +B1*B3  +B2*B3  -B4SQ  -B5SQ  -B6SQ
      A2          = -B1  -B2  -B3
!
!     COMPUTE THE THREE EIGENVALUES GIVEN BY THE ROOTS OF THE 
!     CUBIC POLYNOMIAL. THESE EIGENVALUES SHOULD ALL BE REAL AND 
!     POSITIVE. THE SINGULAR VALUES ARE GIVEN BY THE SQUARE ROOTS 
!     OF THE EIGENVALUES.
!
      SLIDE        = A2*THIRD
      SLIDSQ       = SLIDE*SLIDE
      SLIDCB       = SLIDSQ*SLIDE
      Q            = A1*THIRD  -SLIDSQ
      R            = .5*THIRD*(A1*A2  -3.*A0)  -SLIDCB
      DSQ          = -Q**3  -R*R
      IF (DSQ+TOL.LT.0.) WRITE (6,900) A0,A1,A2,R,Q,DSQ
      DSQ          = ABS(DSQ)
      IF (DSQ.LT.TOL) THEN
         IF (R.LT.TOL) THEN
            S             = 0.
         ELSE
            S             = R**THIRD
         ENDIF
         SIG1(L)       = 2.*S  -SLIDE
         SIG2(L)       = -S  -SLIDE
         SIG3(L)       = SIG2(L)
      ELSE
         D             = SQRT(DSQ)
         THETA         = ATAN2(D,R)
         CT            = COS(THETA*THIRD)
         ST            = SIN(THETA*THIRD)
         QA            = MAX(0.0D0,-Q)
         QROOT         = SQRT(QA)
         P1            = QROOT*CT
         P2            = QROOT*ST*SQRT(3.0D0)
         SIG1(L)       = 2.*P1  -SLIDE
         SIG2(L)       = -P1  -P2  -SLIDE
         SIG3(L)       = -P1  +P2  -SLIDE
      ENDIF
      IF (MIN(SIG1(L),SIG2(L),SIG3(L)).LT.0.) THEN
         IF (MIN(SIG1(L),SIG2(L),SIG3(L))+TOL.LT.0.) THEN
            GO TO 320
         ELSE
            SIG1(L)     = MAX(0.0D0,SIG1(L))
            SIG2(L)     = MAX(0.0D0,SIG2(L))
            SIG3(L)     = MAX(0.0D0,SIG3(L))
         ENDIF
      ENDIF
      SIG1(L)      = MAX(TOL,SIG1(L))
      SIG2(L)      = MAX(TOL,SIG2(L))
      SIG3(L)      = MAX(TOL,SIG3(L))
      SIG1(L)      = SQRT(SIG1(L))
      SIG2(L)      = SQRT(SIG2(L))
      SIG3(L)      = SQRT(SIG3(L))
      SIGLRG       = MAX(SIG1(L),SIG2(L),SIG3(L))
      SIGLOW       = MIN(SIG1(L),SIG2(L),SIG3(L))
      IF (SIGLOW.GT.1.E-8) THEN
         COND          = SIGLRG/SIGLOW
      ELSE
         COND      = 0.
      ENDIF
      IF (NPASS.EQ.0) THEN
         NPASS         = 1
         CNDMIN        = COND
         CNDMAX        = COND
         SIGMN         = SIGLOW
         SIGMX         = SIGLRG
         LMIN          = L
         LMAX          = L
      ELSE
         IF (SIGLOW.LT.SIGMN) LMIN = L
         IF (SIGLRG.GT.SIGMX) LMAX = L
         SIGMN         = MIN(SIGMN,SIGLOW)
         SIGMX         = MAX(SIGMX,SIGLRG)
         CNDMIN        = MIN(CNDMIN,COND)
         CNDMAX        = MAX(CNDMAX,COND)
      ENDIF
  100 CONTINUE
      RETURN
  300 WRITE (6,600) N1,N2,N3,N4
      STOP
  310 WRITE (6,610) 
      WRITE (6,900) A0,A1,A2,R,Q,DSQ
      STOP
  320 WRITE (6,620) L,SIG1(L),SIG2(L),SIG3(L)
      STOP
  600 FORMAT(//5X,'A ZERO VOLUME TETRAHEDRON HAS BEEN FOUND IN THE '// &
                  'UNDEFORMED MESH'//5X,'THE VERTICES OF THE ELEMENT '// &
                  'ARE ',I7,' , ',I7,' , ',I7,' AND ',I7// &
               5X,'PROGRAM STOPPED IN ROUTINE SNGLAR.')
  610 FORMAT(//5X,'THE DISCRIMINANT ASSOCIATED WITH THE SOLUTION OF', &
               5X,'CUBIC EQUATION FOR THE SINGULAR VALUES OF THE', &
               5X,'DEFORMATION MATRIX IS POSITIVE, THUS INDICATING', &
               5X,'ONLY ONE REAL ROOT. THERE APPEARS TO BE AN ERROR IN', &
               5X,'THE CODING FOR THE COMPUTATION OF SINGULAR VALUES.', &
               5X,'PROGRAM STOPPED IN ROUTINE SNGLAR.')
  620 FORMAT(//5X,'AT LEAST ONE OF THE SINGULAR VALUES SQUARED IS', &
                5X,'NEGATIVE. THERE THUS APPEARS TO BE AN ERROR IN', &
               5X,'THE CODING FOR THE COMPUTATION OF SINGULAR VALUES.', &
               5X,'PROGRAM STOPPED IN ROUTINE SNGLAR.'//5X,'CELL ',I7, &
               ' SIG1 = ',E13.5,' SIG2 = ',E13.5,' SIG3 = ',E13.5)
  900 FORMAT(/5X,'A POSITIVE DISCRIMINANT HAS BEEN COMPUTED WHEN', &
              5X,'SOLVING THE CUBIC EQUATION FOR THE SINGULAR VALUES', &
              5X,'OF THE DEFORMATION MATRIX, THUS INDICATING ONLY', &
              5X,'ONE REAL ROOT.'// 5X,'COEFFICIENTS ARE ', &
                 'A0 = ',F10.4,' A1 = ',F10.4,' A2 = ',F10.4, &
              5X,'R = ',E13.5,' Q = ',E13.5,' DSQ = ',E13.5)
      END SUBROUTINE SNGLAR




!
!     ******************************************************************
!
      SUBROUTINE DENSFN (X,NNODE,NFCE,NBFACE,NEDGE,NDG,DENS,RESID,FVCNT)
!
!     ******************************************************************
!     *                                                                *
!     *  EVALUATE DENSITY FUNCTION FOR ALL MESH POINTS                 *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1998                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: NBFACE,NEDGE,NNODE
      INTEGER :: NDG(2,*),NFCE(3,*)
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: FVCNT(*),RESID(*)
      
      INTEGER :: I,IT,J,K,N,N1,N2
      DOUBLE PRECISION :: DIFDEN,DIST,EPS,RMAX,RMAX1
!
!     ******************************************************************
!
      EPS       = 1.
!
!     SET INITIAL DENSITY FUNCTION TO ZERO
!
      DO 10 N=1,NNODE
      FVCNT(N)  = 0.
      DENS(N)   = 0.
   10 CONTINUE
!
!     INITIALIZE FVCNT TO -1 FOR BOUNDARY POINTS
!
      DO 15 J=1,NBFACE
      DO 15 K=1,3
      N         = NFCE(K,J)
      IF (N.GT.0) THEN
         FVCNT(N)  = -1.
      ENDIF
   15 CONTINUE
!
!     ESTABLISH INITIAL GUESS FOR DENSITY FUNCTION FOR ALL MESH POINTS
!
      DO 20 I=1,NEDGE
      N1        = NDG(1,I)
      N2        = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 20
      DIST      = SQRT((X(1,N1)  -X(1,N2))**2  +(X(2,N1)  -X(2,N2))**2 &
                     +(X(3,N1)  -X(3,N2))**2)
      IF (FVCNT(N1).GE.0.) THEN
         FVCNT(N1) = FVCNT(N1)  +1.
         DENS(N1)  = DENS(N1)  +DIST
      ELSE
         IF (FVCNT(N2).LT.0.) THEN
            FVCNT(N1) = FVCNT(N1)  -1.
            DENS(N1)  = DENS(N1)  +DIST
            FVCNT(N2) = FVCNT(N2)  -1.
            DENS(N2)  = DENS(N2)  +DIST
         ENDIF
      ENDIF
      IF (FVCNT(N2).GE.0.) THEN
         FVCNT(N2) = FVCNT(N2)  +1.
         DENS(N2)  = DENS(N2)  +DIST
      ENDIF
   20 CONTINUE
!
!     NORMALIZE DENSITY FUNCTION BY NUMBER OF EDGES MEETING AT A VERTEX
!
      DO 30 N=1,NNODE
      IF (FVCNT(N).LT.0.) FVCNT(N) = FVCNT(N)  +1.
      IF (ABS(FVCNT(N)).LT.1.) GO TO 30
      FVCNT(N)  = 1./FVCNT(N)
      DENS(N)   = ABS(FVCNT(N))*DENS(N)
   30 CONTINUE
!
!     START OF ITERATIVE CYCLE
!
      IT        = 0
   65 IT        = IT  +1
!
!     SET INITIAL RESIDUALS TO ZERO
!
      DO 70 N=1,NNODE
      RESID(N)  = 0.
   70 CONTINUE
!
!     ACCUMMULATE EDGE DIFFERENCES OF DENSITY FUNCTION FOR ALL
!     EDGES INCIDENT TO EACH POINT
!
      DO 75 I=1,NEDGE
      N1        = NDG(1,I)
      N2        = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 75
      DIFDEN    = DENS(N1)  -DENS(N2)
      RESID(N1) = RESID(N1)  -DIFDEN
      RESID(N2) = RESID(N2)  +DIFDEN
   75 CONTINUE
!
!     UPDATE DENSITY FUNCTION AT EACH NON-BOUNDARY POINT 
!
      RMAX      = 0.
      DO 80 N=1,NNODE
      IF (FVCNT(N).LT.0.) GO TO 80
      DENS(N)   = DENS(N)  +EPS*FVCNT(N)*RESID(N)
      RMAX      = MAX(RMAX,ABS(RESID(N)))
   80 CONTINUE

!     WRITE (6,900) IT,RMAX

      IF (IT.EQ.1) RMAX1 = RMAX
      IF (RMAX.GT.MAX(0.0001*RMAX1,1.0D0-9).AND.IT.LT.1000) GO TO 65
      WRITE (6,900) IT,RMAX
  900 FORMAT('IN DENSFN...,  IT = ',I4,' RMAX = ',F10.6)
!
!     ITERATIVE CYCLE IS COMPLETE
!
!      TIM       = SECOND (0)
!      WRITE (6,700) TIM
!  700 FORMAT(/'TIME = ',F7.1,' SECONDS'/)
      RETURN
      END SUBROUTINE DENSFN





!
!     ******************************************************************
!
      SUBROUTINE EDGLEN (X,NNODE,ITYP,NEDGE,NDG,DENS,NVCNT)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPARE EDGE LENGTHS WITH DENSITY FUNCTION.                   *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   2001                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: NEDGE,NNODE
      INTEGER :: NDG(2,*),ITYP(*),NVCNT(*)
      DOUBLE PRECISION :: X(3,*),DENS(*)
      
      INTEGER :: I,N,NHIGH,NLOW,NMID,NPASS,NTOT,NVMAX,NVMIN,N1,N2
      DOUBLE PRECISION :: DENSAV,DSQAV
!
!     ******************************************************************
!
      NLOW      = 0
      NMID      = 0
      NHIGH     = 0
!
!     COMPUTE EDGE VALENCE OF EACH MESH POINT
!
      DO 10 N=1,NNODE
      NVCNT(N)  = 0
   10 CONTINUE
      DO 20 I=1,NEDGE
      N1        = NDG(1,I)
      N2        = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 20
      NVCNT(N1) = NVCNT(N1)  +1
      NVCNT(N2) = NVCNT(N2)  +1
      DENSAV    = .5*(DENS(N1)  +DENS(N2))
      DSQAV     = SQRT((X(1,N1)  -X(1,N2))**2  +(X(2,N1)  -X(2,N2))**2 &
                      +(X(3,N1)  -X(3,N2))**2)
      IF (DSQAV.LT.0.5*DENSAV) NLOW  = NLOW  +1
      IF (DSQAV.GT.1.5*DENSAV) NHIGH = NHIGH  +1
      IF (DSQAV.GE.0.5*DENSAV.AND.DSQAV.LE.1.5*DENSAV) NMID = NMID  +1
   20 CONTINUE
      NTOT      = NLOW  +NMID  +NHIGH
      WRITE (6,600) NLOW,NMID,NHIGH,NTOT
      NVMAX     = 0
      NVMIN     = 0
      NPASS   = 0
      DO 25 N=1,NNODE
      IF (ITYP(N).GE.0) THEN
         IF (NPASS.EQ.0) THEN
            NPASS   = 1
            NVMAX   = NVCNT(N)
            NVMIN   = NVCNT(N)
         ELSE
            NVMAX   = MAX(NVMAX,NVCNT(N))
            NVMIN   = MIN(NVMIN,NVCNT(N))
         ENDIF
!     ELSE
!        WRITE (6,920) N,ITYP(N)
! 920    FORMAT('pt ',I5,' ITYP ',I3)
      ENDIF
   25 CONTINUE
      WRITE (6,610) NVMIN,NVMAX
      RETURN
  600 FORMAT(/5X,'NLOW = ',I6,' NMID = ',I6,' NHIGH = ',I6, &
                 ' NTOT = ',I6)
  610 FORMAT(/5X,'MIN EDGE VALENCE = ',I3,' MAX EDGE VALENCE = ',I3)
      END SUBROUTINE EDGLEN





!
!     ******************************************************************
!
      SUBROUTINE COARSN (LC,X,NNODE,NDC,NBH,IPROT,NCELL, &
                         ITYP,XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                         NVCNT,DENS,IFLAG,NFLAG,NPTET, &
                         NEDGE,NDG,IDGP,NDGP, &
                         NOCTR,IOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE,NREF, &
                         FAC,SIG1,SIG2,SIG3, &
                         KSRCH,NSRCH,IRING,NTETKP, &
                         LNBR,ISHK,MNBR,KSHK,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  CALL ROUTINE COLAPS TO MERGE POINTS ASSOCIATED WITH SELECTED  *
!     *  EDGES                                                         *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1998                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,LC,NCELL,NEDGE,NNODE
      INTEGER :: ITYP(*),NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*)
      INTEGER :: IDGP(*),NDGP(*),NVCNT(*),NPTET(*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: IFLAG(*),NFLAG(*),NSRCH(*),KSRCH(*)
      INTEGER :: IRING(*),NTETKP(*),LNBR(*),ISHK(*),MNBR(*),KSHK(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*),FAC(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: SIG1(*),SIG2(*),SIG3(*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)
      
      INTEGER :: I,J,L,N,NCOL,NFAIL,NPTS,NTEST,NTET,N1,N2,N3,N4
      DOUBLE PRECISION :: DENSAV,DIST,SIGMAX,SIGMIN,EXCESS,EXPAND
!
!     ******************************************************************
!
!
!     COMPUTE EDGE VALENCE OF EACH MESH POINT
!
      DO 10 N=1,NNODE
      NVCNT(N)  = 0
      FAC(N)    = 0.
   10 CONTINUE
      DO 20 I=1,NEDGE
      N1        = NDG(1,I)
      N2        = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 20
      NVCNT(N1) = NVCNT(N1)  +1
      NVCNT(N2) = NVCNT(N2)  +1
   20 CONTINUE
!
!     SELECT TETRAHEDRA TO BE REMOVED
!
      NCOL    = 0
      DO 30 J=1,NCELL
      IF (NBH(1,J).LE.0) GO TO 30
      IF (NDC(4,J).EQ.-1) GO TO 30
      IF (LC.EQ.1) THEN
         SIGMAX  = MAX(SIG1(J),SIG2(J),SIG3(J))
         SIGMIN  = MIN(SIG1(J),SIG2(J),SIG3(J))
         EXPAND  = MIN(1.0D0,SIGMAX/SIGMIN  -1.)
      ELSE
         EXCESS  = MAX(0.0D0,RAT(J)-20.)
         EXPAND  = MIN(1.0D0,EXCESS/100.)
      ENDIF
      N1      = NDC(1,J)
      N2      = NDC(2,J)
      N3      = NDC(3,J)
      N4      = NDC(4,J)
      FAC(N1) = FAC(N1)  +EXPAND
      FAC(N2) = FAC(N2)  +EXPAND
      FAC(N3) = FAC(N3)  +EXPAND
      FAC(N4) = FAC(N4)  +EXPAND
   30 CONTINUE
      DO 35 N=1,NNODE
      IF (ITYP(N).GE.0) THEN
         DENS(N) = DENS(N)*(1.  +FAC(N)/FLOAT(NVCNT(N)))
      ENDIF
   35 CONTINUE
!
!     COLLAPSE EDGES THAT ARE TOO SHORT COMPARED WITH THE LENGTH
!     DENSITY FUNCTION.
!
      DO 40 I=1,NEDGE
      N1        = NDG(1,I)
      N2        = NDG(2,I)
      IF (N1.LE.0.OR.N2.LE.0) GO TO 40
      DENSAV    = .5*(DENS(N1)  +DENS(N2))
      DIST      = SQRT((X(1,N1)  -X(1,N2))**2  +(X(2,N1)  -X(2,N2))**2 &
                      +(X(3,N1)  -X(3,N2))**2)
      IF (DIST.GE.0.4*DENSAV) GO TO 40

      CALL COLAPS (N1,N2,NCOL,NTEST,NFAIL, &
                  X,NDC,NBH,IPROT, &
                  ITYP,XCEN,YCEN,ZCEN,VOL,RC,RAT, & 
                  NVCNT,IFLAG,NFLAG,NPTET, &
                  NDG,IDGP,NDGP, &
                  NOCTR,IOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE,NREF, &
                  KSRCH,NSRCH,IRING,NTETKP, &
                  LNBR,ISHK,MNBR,KSHK,TOLV)

   40 CONTINUE
!
!     WRITE OUT NUMBER OF MESH POINTS AND TRIANGLES
!
      NPTS     = 0
      DO 50 N=1,NNODE
      IF (ITYP(N).LT.0) GO TO 50
      NPTS     = NPTS  +1
   50 CONTINUE
      NTET     = 0
      DO 60 L=1,NCELL
      IF (NBH(1,L).EQ.0) GO TO 60
      IF (NDC(4,L).EQ.-1) GO TO 60
      NTET     = NTET  +1
   60 CONTINUE
      WRITE (6,730) NPTS,NCOL,NTET
      RETURN
  730 FORMAT(/'MESH COARSENING COMPLETE'/ &
           5X,I7,' TOTAL MESH POINTS'/ &
           5X,I7,' FIELD POINTS REMOVED'/ &
           5X,I7,' MESH CELLS'/)
      END SUBROUTINE COARSN







!
!     ******************************************************************
!
      SUBROUTINE RADRAT (X,NNODE,NDC,NBH,IPROT,NCELL,NFCE,NBFACE, &
                         ITYP,IPOINT,VOL,RC,RAT)
!
!     ******************************************************************
!     *                                                                *
!     *  DETERMINE MINIMUM AND MAXIMUM VOLUMES AND STATISTICS ON THE   *
!     *  RATIO OF CIRCUMRADIUS TO IN-RADIUS FOR ALL TETRAHEDRA         *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!  
      IMPLICIT NONE

      INTEGER :: NBFACE,NCELL,NNODE
      INTEGER :: NBH(4,*),IPROT(*),ITYP(*),IPOINT(*)
      INTEGER :: NFCE(3,*),NDC(4,*)
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: VOL(*),RC(*),RAT(*)
!
      INTEGER :: J,JAMIN,JMAX,JMIN,K,L,N,NBFC,NCNT,NSURPT,NVOLPT, &
                 N1,N2,N3
      DOUBLE PRECISION :: V(4),AR(4)
      DOUBLE PRECISION :: ANGL,ANGL1,ANGL2,ANGL3,ANGMAX,ANGMIN,DISMIN, &
                          FAREA,FMAX,FMIN,Q,QMAX,QMIN,VMAX,VMIN, &
                          RATAVR,RATMAX,RATMIN,RCMAX,RCMIN,SIGMA,SIGSQ
!
!     ******************************************************************
!
      RATAVR    = 0.
      NCNT      = 0
      DO 20 L=1,NCELL
      IF (IPROT(L).EQ.1) GO TO 20
      IF (NBH(1,L).EQ.0) GO TO 20
      IF (NDC(4,L).LE.0) GO TO 20
      NCNT      = NCNT  +1
      RATAVR    = RATAVR  +RAT(L)
      IF (NCNT.GT.1) GO TO 10
      VMIN      = VOL(L)
      VMAX      = VOL(L)
      RCMIN     = RC(L)
      RCMAX     = RC(L)
      RATMIN    = RAT(L)
      RATMAX    = RAT(L)
      GO TO 20
   10 VMIN      = MIN(VOL(L),VMIN)
      VMAX      = MAX(VOL(L),VMAX)
      RCMIN     = MIN(RCMIN,RC(L))
      RCMAX     = MAX(RCMAX,RC(L))
      RATMIN    = MIN(RATMIN,RAT(L))
      RATMAX    = MAX(RATMAX,RAT(L))
   20 CONTINUE
      RATAVR    = RATAVR/FLOAT(NCNT)
      VMIN      = VMIN/6.0D0
      VMAX      = VMAX/6.0D0
!
!     COMPUTE STANDARD DEVIATION OF CIRCUMRADIUS TO IN-RADIUS
!     RATIO
!
      SIGSQ     = 0.
      DO 30 L=1,NCELL
      IF (IPROT(L).EQ.1) GO TO 30
      IF (NBH(1,L).EQ.0) GO TO 30
      IF (NDC(4,L).LE.0) GO TO 30
      SIGSQ     = SIGSQ  +(RAT(L)  -RATAVR)**2
   30 CONTINUE
      SIGMA     = SQRT(SIGSQ)/FLOAT(NCNT)
      WRITE (6,600) VMIN,VMAX,RCMIN,RCMAX,RATMIN,RATMAX,RATAVR,SIGMA

      CALL TETANG (X,NDC,NBH,IPROT,NCELL)
!
!     COMPUTE AREA AND ANGLES OF EACH SURFACE TRIANGLE
!
      JMAX        = 1
      JMIN        = 1
      JAMIN       = 1
      ANGMIN      = 360.
      ANGMAX      = 0.
      QMIN        = 1.E15
      QMAX        = 1.
      FMAX        = -1.
      FMIN        = -1.
      NBFC        = 0
      DO 75 J=1,NBFACE
      N1          = NFCE(1,J)
      N2          = NFCE(2,J)
      N3          = NFCE(3,J)
      IF (N1.LT.0) GO TO 75
      NBFC        = NBFC  +1

      CALL FANGLE (J,X,NFCE,ANGL1,ANGL2,ANGL3,Q)

      ANGL        = MIN(ANGL1,ANGL2,ANGL3)
      IF (ANGL.GT.ANGMIN) GO TO 60
      JAMIN       = J
      ANGMIN      = ANGL
   60 ANGMAX      = MAX(ANGL1,ANGL2,ANGL3,ANGMAX)
      QMIN        = MIN(Q,QMIN)
      QMAX        = MAX(Q,QMAX)
      FAREA       = FACEAR (X,N1,N2,N3)
      IF (FMAX.LT.0.) FMAX = FAREA
      IF (FMIN.LT.0.) FMIN = FAREA
      IF (FAREA.LT.FMAX) GO TO 70
      JMAX        = J
      FMAX        = FAREA
   70 IF (FAREA.GT.FMIN) GO TO 75
      JMIN        = J
      FMIN        = FAREA
   75 CONTINUE
      WRITE (6,610) ANGMIN,ANGMAX,QMIN,QMAX
!
!     COUNT NUMBER OF POINTS IN THE MESH AND ON THE BOUNDARY SURFACE
!
      NVOLPT      = 0
      NSURPT      = 0
      DO 80 N=1,NNODE
      IPOINT(N)   = 0
      IF (ITYP(N).GE.0) NVOLPT = NVOLPT  +1
   80 CONTINUE
      DO 85 J=1,NBFACE
      DO 85 K=1,3
      N           = NFCE(K,J)
      IF (IPOINT(N).EQ.0) THEN
         IPOINT(N)   = 1
         NSURPT      = NSURPT  +1
      ENDIF
   85 CONTINUE
      DO 90 N=1,NNODE
      IPOINT(N)   = 0
   90 CONTINUE
      WRITE (6,620) NSURPT,NBFC,NVOLPT,NCNT
      RETURN
  600 FORMAT(//5X,'*************************************************'/ &
               5X,'*                                               *'/ &
               5X,'*  STATISTICS OF MESH TETRAHEDRA                *'/ &
               5X,'*                                               *'/ &
               5X,'*  MINIMUM VOLUME = ',E13.5,'               *'/ &
               5X,'*  MAXIMUM VOLUME = ',E13.5,'               *'/ &
               5X,'*                                               *'/ &
               5X,'*  MINIMUM CIRCUMRADIUS = ',E13.5,'         *'/ &
               5X,'*  MAXIMUM CIRCUMRADIUS = ',E13.5,'         *'/ &
               5X,'*                                               *'/ &
               5X,'*  MINIMUM RADIUS RATIO = ',E13.5,'         *'/ &
               5X,'*  MAXIMUM RADIUS RATIO = ',E13.5,'         *'/ &
               5X,'*  AVERAGE RADIUS RATIO = ',E13.5,'         *'/ &
               5X,'*  STANDARD DEVIATION   = ',E13.5,'         *'/ &
               5X,'*                                               *'/ &
               5X,'*************************************************')
  610 FORMAT(  5X,'*************************************************'/ &
               5X,'*                                               *'/ &
               5X,'*  STATISTICS OF BOUNDARY SURFACE TRIANGLES     *'/ &
               5X,'*                                               *'/ &
               5X,'*  MINIMUM BOUNDARY FACE ANGLE = ',F6.2,'         *'/ &
               5X,'*  MAXIMUM BOUNDARY FACE ANGLE = ',F6.2,'         *'/ &
               5X,'*  MINIMUM BOUNDARY RADIUS RATIO = ',F6.2,'       *'/ &
               5X,'*  MAXIMUM BOUNDARY RADIUS RATIO = ',F6.2,'       *'/ &
               5X,'*                                               *'/ &
               5X,'*************************************************')
  620 FORMAT(  5X,'*************************************************'/ &
               5X,'*                                               *'/ &
               5X,'*  TOTAL NUMBER OF SURFACE POINTS = ',I6,'      *'/ &
               5X,'*  TOTAL NUMBER OF SURFACE FACES  = ',I6,'      *'/ &
               5X,'*  TOTAL NUMBER OF MESH POINTS    = ',I6,'      *'/ &
               5X,'*  TOTAL NUMBER OF TETRAHEDRA     = ',I6,'      *'/ &
               5X,'*                                               *'/ &
               5X,'*************************************************')
      END SUBROUTINE RADRAT






!
!     ******************************************************************
!
      SUBROUTINE VOLPUT (X,ITYP,NBPTS,NNODE,NDC,NBH,IPROT,NCELL, &
                         NDG,IDGP,NDGP,NEDGE, &
                         VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                         IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                         XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                         XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                         IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                         NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                         NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                         JLAST,JFIRST,NTRACK,VOLMIN,RCMX,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  INSERTION OF FIELD POINTS INTO DELAUNAY TRIANGULATION         *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,JFIRST,JLAST,NBPTS,NCELL,NEDGE,NNODE,NTRACK
      INTEGER :: NDC(4,*),NDG(2,*),NBH(4,*),IPROT(*),IDGP(*),NDGP(*)
      INTEGER :: ITYP(*),NPTET(*),NACPT(*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: NPOINT(*),IPOINT(*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: LNKUP(*),LNKDN(*)
      INTEGER :: NSRCH(*),KSRCH(*)
      INTEGER :: NTRI(3,*),NFILL(*),NEWNBH(4,*),NOLD(*), &
                 NEWCEL(*),NSHAKE(*),NCAV(4,*)
      INTEGER :: NEDGRM(*)
      INTEGER :: LDEL(*),NCAVFC(3,*),IKEEP(*)      
      DOUBLE PRECISION :: VOLMIN,RCMX,TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XOCTR(2,*),YOCTR(2,*),ZOCTR(2,*), &
                          XHOLD(2,*),YHOLD(2,*),ZHOLD(2,*), &
                          XKEEP(2),YKEEP(2),ZKEEP(2),XFAR(2), & 
                          YFAR(2),ZFAR(2)
      DOUBLE PRECISION :: XC(*),YC(*),ZC(*),V(*),RAD(*),RCRIN(*)

      INTEGER :: J,JPAST,JPRE,K,K1,K2,L,LBRK,LNBH,N,NCLOSE,NCYC,NDIFF, &
                 NFAIL,NINP,NPASS,NPTS,NSRFPT,NSTART,NTET,N1,N2,N3,N4,N5
      DOUBLE PRECISION :: DENT,DISMIN,RNX,RNY,RNZ,SCLDIS,VTET1,VTET2, &
                          XPT,YPT,ZPT
!
!     ******************************************************************
!
!     SET INTERIOR DISTANCE SCALING
!
      SCLDIS    = .75
!
!     INITIALIZE INTEGER VARIABLES
!
      NSRFPT    = NBPTS  +8
      NSTART    = NNODE
      NCYC      = 0
      NINP      = 0
!
!     SET UP LINKED LIST OF CELLS
!
      JLAST    = 0
      JFIRST   = 0
      NTRACK   = 0
      J        = 0
   10 J        = J  +1
      IF (J.GT.NCELL) GO TO 15
      IF (IPROT(J).EQ.1) GO TO 10
      IF (NBH(1,J).LE.0) GO TO 10
      NACPT(J) = 0
      N1       = NDC(1,J)
      N2       = NDC(2,J)
      N3       = NDC(3,J)
      N4       = NDC(4,J)
      DENT     = DENS(N1)  +DENS(N2)  +DENS(N3)  +DENS(N4)
      IF (RC(J).LT.0.22*DENT) GO TO 10
      NACPT(J) = 1
      NTRACK   = NTRACK  +1
      IF (JLAST.NE.0) LNKDN(JLAST) = J
      LNKDN(J) = 0
      LNKUP(J) = JLAST
      JLAST    = J
      IF (JFIRST.NE.0) GO TO 10
      JFIRST   = J
      GO TO 10
   15 WRITE (6,910) NTRACK
  910 FORMAT('IN VOLPUT,  NTRACK = ',I6)
      IF (NTRACK.EQ.0) GO TO 70
!
!     START OF ITERATIVE LOOP FOR POINT INSERTION
!
   30 NCYC      = NCYC  +1
      J         = JFIRST
      IF (NACPT(J).EQ.0) GO TO 215
      IF (IPROT(J).EQ.1) GO TO 200
      XPT       = XCEN(J)
      YPT       = YCEN(J)
      ZPT       = ZCEN(J)
      IF (XPT.LT.XFAR(1).OR.XPT.GT.XFAR(2)) GO TO 45
      IF (YPT.LT.YFAR(1).OR.YPT.GT.YFAR(2)) GO TO 45 
      IF (ZPT.LT.ZFAR(1).OR.ZPT.GT.ZFAR(2)) GO TO 45
      DO 32 K=1,4
      LNBH      = NBH(K,J)
      IF (IPROT(LNBH).EQ.0) GO TO 32

      CALL LOCK (J,LNBH,N1,N2,N3,N4,N5,K1,K2,NDC,NBH)

      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      VTET1   = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                         +RNZ*(X(3,N4)  -X(3,N1))
      VTET2   = RNX*(XPT  -X(1,N1))  +RNY*(YPT  -X(2,N1)) &
                                     +RNZ*(ZPT  -X(3,N1))
      IF (VTET1*VTET2.LT.0.) GO TO 45
   32 CONTINUE
      NNODE     = NNODE  +1
      IF (NNODE.GT.MXNODE) GO TO 230
      NPASS     = 0
      N         = NNODE
      X(1,N)    = XPT
      X(2,N)    = YPT
      X(3,N)    = ZPT
   34 ITYP(N)   = 8
      IDONE(N)  = 0
      IPOINT(N) = 0

      CALL OCTFND (0,NCLOSE,X(1,N),X(2,N),X(3,N),0,0,0,DISMIN, &
                   X,NOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE, &
                   XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                   XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH)

      IF (DISMIN.LT.1.E-15) GO TO 40
      IF (NCLOSE.LE.NSRFPT) GO TO 35
      IF (NCLOSE.LE.NSTART) GO TO 35
      IF (SQRT(DISMIN).LT.SCLDIS*DENS(NCLOSE)) GO TO 40
   35 IDONE(N)  = 1
      LBRK      = NPTET(NCLOSE)

      CALL INSERT (N,NCLOSE,LBRK,J,NFAIL, &
                   X,NDC,NBH,IPROT,NCELL,NDG,IDGP,NDGP,NEDGE, &
                   VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                   IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                   NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                   NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                   JLAST,JFIRST,NTRACK,VOLMIN,RCMX,TOLV)

      IF (NFAIL.EQ.5.AND.NPASS.EQ.0) THEN
         N1         = NDC(1,J)
         N2         = NDC(2,J)
         N3         = NDC(3,J)
         N4         = NDC(4,J)
         XPT        = .25*(X(1,N1)  +X(1,N2)  +X(1,N3)  +X(1,N4))
         YPT        = .25*(X(2,N1)  +X(2,N2)  +X(2,N3)  +X(2,N4))
         ZPT        = .25*(X(3,N1)  +X(3,N2)  +X(3,N3)  +X(3,N4))
         NPASS      = 1
         GO TO 34
      ENDIF
      IF (NFAIL.GE.1) GO TO 40
!
!     POINT ACCEPTED AND INSERTED INTO MESH
!

!C     WRITE (6,880) N,X(1,N),X(2,N),X(3,N),DENS(N)
!C 880 FORMAT('+++ N ',I6,' X,Y,Z ',3(2X,F8.3),' DENS ',F8.3)

      NINP      = NINP  +1

      CALL OCTFIL (N,X,NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR)

      GO TO 50
!
!     POINT NOT ACCEPTED FOR INSERTION
!
   40 NNODE     = NNODE  -1
   45 NACPT(J)  = 0
      JPRE      = LNKUP(J)
      JPAST     = LNKDN(J)
      IF (JPRE.NE.0) LNKDN(JPRE) = JPAST
      IF (JPAST.NE.0) LNKUP(JPAST) = JPRE
      IF (JPRE.EQ.0) JFIRST = JPAST
      IF (JPAST.EQ.0) JLAST = JPRE
      NTRACK    = NTRACK  -1
   50 IF (MOD(NCYC,10000).EQ.0) WRITE (6,710) NCYC,NINP,NTRACK
      IF (JFIRST.NE.0) GO TO 30
!
!     ITERATIVE INSERTION OF POINTS IS COMPLETE
!
   70 WRITE (6,710) NCYC,NINP,NTRACK
      NPTS     = 0
      DO 75 N=1,NNODE
      IF (N.GT.0) THEN
         IF (ITYP(N).GE.0) NPTS = NPTS  +1
      ENDIF
   75 CONTINUE
      NDIFF    = NNODE  -NSTART
      NTET     = 0
      DO 80 L=1,NCELL
      IF (IPROT(L).EQ.1) GO TO 80
      IF (NBH(1,L).EQ.0) GO TO 80
      NTET     = NTET  +1
   80 CONTINUE
      WRITE (6,720) NPTS,NDIFF,NTET
      RETURN
  200 WRITE (6,600)
      STOP
  210 WRITE (6,610) NFAIL,N,X(1,N),X(2,N),X(3,N)
      STOP
  215 WRITE (6,615)
      STOP
  230 WRITE (6,630)
      STOP
  710 FORMAT(/' NCYC = ',I7,' POINTS INSERTED = ',I7,' NTRACK = ',I7)
  720 FORMAT(/'ADAPTIVE REFINEMENT COMPLETE'/ &
            5X,I7,' TOTAL MESH POINTS'/ &
            5X,I7,' FIELD POINTS INSERTED'/ &
            5X,I7,' MESH CELLS'/)
  600 FORMAT(//'A PROTECTED TETRAHEDRON IS AMONG THOSE TO BE REFINED')
  610 FORMAT(//5X,'NFAIL = ',I6,' N ',I6,' X,Y,Z ',3(2X,F6.2)/ &
               5X,'PROGRAM STOPPED IN VOLPUT')
  615 FORMAT(//5X,'LABEL NACPT IS ZERO FOR CELL FROM ACTIVE LIST'/ &
               5X,'PROGRAM STOPPED IN VOLPUT')
  630 FORMAT(//'NUMBER OF POINTS INSERTED EXCEEDS DIMENSION OF ARRAY X'/ &
               'INCREASE SIZE OF MNODE.'/ &
               'PROGRAM STOPPED IN ROUTINE VOLPUT.')
      END SUBROUTINE VOLPUT





!
!     ******************************************************************
!
      SUBROUTINE PUTPNT (J,NFAIL, &
                         X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                         NDG,IDGP,NDGP,NEDGE, &
                         VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                         IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                         XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                         XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                         IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                         NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                         NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                         VOLMIN,RCMX,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  INSERTION OF A POINT INTO DELAUNAY TRIANGULATION TO ELIMINATE *
!     *  A SLIVER LIKE TETRAHEDRON.                                    *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   2001                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,J,NCELL,NEDGE,NFAIL,NNODE
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: NTRI(3,*),NFILL(*),NEWNBH(4,*),NOLD(*),NEWCEL(*), &
                 NSHAKE(*),NCAV(4,*)
      INTEGER :: NEDGRM(*),LDEL(*),NCAVFC(3,*),IKEEP(*)
      INTEGER :: ITYP(*),NPTET(*),NACPT(*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: NPOINT(*),IPOINT(*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: LNKUP(*),LNKDN(*)
      INTEGER :: NSRCH(*),KSRCH(*)
      DOUBLE PRECISION :: VOLMIN,RCMX,TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XOCTR(2,*),YOCTR(2,*),ZOCTR(2,*), &
                          XHOLD(2,*),YHOLD(2,*),ZHOLD(2,*), &
                          XFAR(2),YFAR(2),ZFAR(2),XKEEP(2), &
                          YKEEP(2),ZKEEP(2)
      DOUBLE PRECISION :: XC(*),YC(*),ZC(*),V(*),RAD(*),RCRIN(*)

      INTEGER :: JFIRST,JLAST,K,K1,K2,LBRK,LNBH,N,NCLOSE,NPASS,NTRACK, & 
                 N1,N2,N3,N4,N5
      DOUBLE PRECISION :: DISMIN,RNX,RNY,RNZ,VTET1,VTET2,XPT,YPT,ZPT
!
!     ******************************************************************
!
      NFAIL     = 0
      XPT       = XCEN(J)
      YPT       = YCEN(J)
      ZPT       = ZCEN(J)
      IF (XPT.LT.XFAR(1).OR.XPT.GT.XFAR(2)) GO TO 50
      IF (YPT.LT.YFAR(1).OR.YPT.GT.YFAR(2)) GO TO 50 
      IF (ZPT.LT.ZFAR(1).OR.ZPT.GT.ZFAR(2)) GO TO 50
      DO 32 K=1,4
      LNBH      = NBH(K,J)
      IF (IPROT(LNBH).EQ.0) GO TO 32

      CALL LOCK (J,LNBH,N1,N2,N3,N4,N5,K1,K2,NDC,NBH)

      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      VTET1   = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                         +RNZ*(X(3,N4)  -X(3,N1))
      VTET2   = RNX*(XPT  -X(1,N1))  +RNY*(YPT  -X(2,N1)) &
                                     +RNZ*(ZPT  -X(3,N1))
      IF (VTET1*VTET2.LT.0.) GO TO 50
   32 CONTINUE
      NNODE     = NNODE  +1
      IF (NNODE.GT.MXNODE) GO TO 200
      NPASS     = 0
      N         = NNODE
   34 X(1,N)    = XPT
      X(2,N)    = YPT
      X(3,N)    = ZPT
      ITYP(N)   = 8
      IDONE(N)  = 0
      IPOINT(N) = 0

      CALL OCTFND (0,NCLOSE,X(1,N),X(2,N),X(3,N),0,0,0,DISMIN, &
                   X,NOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE, &
                   XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                   XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH)

      IF (DISMIN.LT.1.E-15) GO TO 40
   35 IDONE(N)  = 1
      LBRK      = NPTET(NCLOSE)
      NTRACK    = 0
      JFIRST    = 0
      JLAST     = 0

      CALL INSERT (N,NCLOSE,LBRK,J,NFAIL, &
                   X,NDC,NBH,IPROT,NCELL,NDG,IDGP,NDGP,NEDGE, &
                   VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                   IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                   NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                   NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                   JLAST,JFIRST,NTRACK,VOLMIN,RCMX,TOLV)

      IF (NFAIL.EQ.5.AND.NPASS.EQ.0) THEN
         N1         = NDC(1,J)
         N2         = NDC(2,J)
         N3         = NDC(3,J)
         N4         = NDC(4,J)
         XPT        = .25*(X(1,N1)  +X(1,N2)  +X(1,N3)  +X(1,N4))
         YPT        = .25*(X(2,N1)  +X(2,N2)  +X(2,N3)  +X(2,N4))
         ZPT        = .25*(X(3,N1)  +X(3,N2)  +X(3,N3)  +X(3,N4))
         NPASS      = 1
         GO TO 34
      ENDIF
      IF (NFAIL.GE.1) THEN
         NNODE     = NNODE  -1
         RETURN
      ENDIF
!
!     POINT ACCEPTED AND INSERTED INTO MESH
!
      CALL OCTFIL (N,X,NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR)

      RETURN
!
!     POINT NOT ACCEPTED FOR INSERTION
!
   40 NNODE     = NNODE  -1
   50 NFAIL     = 1
      RETURN
  200 WRITE (6,600)
      STOP
  600 FORMAT(//'NUMBER OF POINTS INSERTED EXCEEDS DIMENSION OF ARRAY X'/ &
               'INCREASE SIZE OF MNODE.'/ &
               'PROGRAM STOPPED IN ROUTINE PUTPNT.')
      END SUBROUTINE PUTPNT






!
!     ******************************************************************
!
      SUBROUTINE INSERT (NP,NCLOSE,LBRK,JCONT,NFAIL, &
                         X,NDC,NBH,IPROT,NCELL,NDG,IDGP,NDGP,NEDGE, &
                         VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                         IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                         NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                         NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                         JLAST,JFIRST,NTRACK,VOLMIN,RCMX,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  POINT (X(1,NP),X(2,NP),X(3,NP)) LIES INSIDE CIRCUMSPHERE OF   *
!     *  TETRAHEDRON LBRK.  FIND THE DELAUNAY CAVITY THAT CONTAINS     *
!     *  THIS POINT AND GENERATE THE DELAUNAY MODIFICATION OF THE      *
!     *  CAVITY TO INCLUDE THE NEW POINT IN THE MESH.                  *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: JCONT,JFIRST,JLAST,LBRK,NCELL,NCLOSE,NEDGE,NFAIL,NP, &
                 NTRACK
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: NPTET(*),NACPT(*)
      INTEGER :: NPOINT(*),IPOINT(*),IFLAG(*),NFLAG(*)
      INTEGER :: NTRI(3,*),NFILL(*),NEWNBH(4,*),NOLD(*), &
                 NEWCEL(*),NSHAKE(*),NCAV(4,*)
      INTEGER :: NEDGRM(*)
      INTEGER :: LNKUP(*),LNKDN(*)
      INTEGER :: LDEL(*),NCAVFC(3,*),IKEEP(*)
      DOUBLE PRECISION :: VOLMIN,RCMX,TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XC(*),YC(*),ZC(*),V(*),RAD(*),RCRIN(*)
      
      INTEGER :: I,IEDGRM,ISMALL,J,JPAST,JPRE,K,L,LCONT,M,MM,NCHK, &
                 NCNT,NCPNT,NDEL,NEDG,NTOT1,NTOT2,N1,N2,N3,N4
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,DENT, &
                          RADMAX,RATMAX,SUM1,SUM2,VDIFF,XPT,YPT,ZPT
!
!     ******************************************************************
!
      NFAIL     = 0
      XPT       = X(1,NP)
      YPT       = X(2,NP)
      ZPT       = X(3,NP)
!
!     FIND TETRAHEDRON LCONT THAT CONTAINS POINT NP
!
      CALL TETLOC (NP,NCLOSE,LBRK,LCONT,NFAIL, &
                   X,NDC,NBH,IPROT,DENS,VOL, &
                  IFLAG,NFLAG,NFILL,NEWCEL,TOLV)

      IF (NFAIL.NE.0) RETURN
!
!     FIND COMPLETE CAVITY OF TETRAHEDRA WHOSE CIRCUMSPHERES CONTAIN
!     THE POINT NP
!
      CALL CAVITY (NP,LCONT,NDEL,LDEL, &
                   X,NDC,NBH,IPROT,IFLAG,NFLAG, &
                   XCEN,YCEN,ZCEN,RC,NFILL,NEWCEL,TOLV)
!
!     DETERMINE FACES ON BOUNDARY OF CAVITY AND GENERATE LIST OF NEW
!     TETRAHEDRA
!
      CALL CAVBND (NP,LCONT,NDEL,LDEL,NCNT, &
                   X,NDC,NBH,VOL,IFLAG,NFLAG, &
                   NFILL,NEWCEL,NTRI,NEWNBH,NOLD,TOLV)
!
!     COMPUTE MAXIMUM CIRCUMRADIUS FOR THE CAVITY TETRAHEDRA
!
      RADMAX    = 0.
      RATMAX    = 0.
      NCHK      = 0
      DO 10 K=1,NDEL
      RADMAX    = MAX(RADMAX,RC(LDEL(K)))
      RATMAX    = MAX(RATMAX,RAT(LDEL(K)))
      IF (LDEL(K).EQ.JCONT) NCHK = 1
   10 CONTINUE
      IF (NCHK.EQ.0) GO TO 295
!
!     COMPUTE VOLUME, CIRCUMCENTER COORDINATES AND CIRCUMRADIUS OF
!     EACH NEW TETRAHEDRON
!
      SUM1      = 0.
      DO 20 K=1,NCNT
      N1        = NP
      N2        = NTRI(1,K)
      N3        = NTRI(2,K)
      N4        = NTRI(3,K)
!
      CALL CIRCUM (X,N1,N2,N3,N4,XC(K),YC(K),ZC(K),V(K),RAD(K), &
                   ISMALL,TOLV)
!
      IF (ISMALL.EQ.1) GO TO 300
      IF (V(K).LT.TOLV) GO TO 310
      IF (RAD(K).GT.RADMAX) GO TO 305
      AREA      = TETAR2 (N1,N2,N3,N4,X)
      RCRIN(K)  = RAD(K)*AREA/V(K)
      IF (RCRIN(K).GT.RATMAX) GO TO 305
      SUM1      = SUM1  +V(K)
   20 CONTINUE
!
!     CHECK PROXIMITY OF POINT TO BOUNDARY SURFACE
!
!!    DO 25 K=1,NCNT
!!    LNBH      = NEWNBH(1,K)
!!    IF (IPROT(LNBH).EQ.0) GO TO 25
!!    N1        = NTRI(1,K)
!!    N2        = NTRI(2,K)
!!    N3        = NTRI(3,K)
!!    RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
!!    RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
!!    RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
!!    FAC       = 1./SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
!!    RNX       = FAC*RNX
!!    RNY       = FAC*RNY
!!    RNZ       = FAC*RNZ
!!    PROJ1     = RNX*(X(1,NP)  -X(1,N1))  +RNY*(X(2,NP)  -X(2,N1))
!!   .                                     +RNZ*(X(3,NP)  -X(3,N1))
!!    PROJ2     = RNX*(XC(K)  -X(1,N1))  +RNY*(YC(K)  -X(2,N1))
!!   .                                   +RNZ*(ZC(K)  -X(3,N1))
!!    PROJ1     = ABS(PROJ1)
!!    IF (6.*PROJ1.LT.(DENS(N1)+DENS(N2)+DENS(N3))) GO TO 320
!! 25 CONTINUE
!
!     CHECK CAVITY VOLUME FOR VISIBILITY
!
      SUM2      = 0.
      DO 30 K=1,NDEL
      SUM2      = SUM2  +VOL(LDEL(K))
   30 CONTINUE
      IF (SUM1.GT.SUM2*(1.+TOLV)) GO TO 350
!
!     UPDATE DATA STRUCTURE FOR FACE AND EDGE LISTS
!
      CALL RECON (LDEL,NDEL,NCNT, &
                  NDC,NBH,IPROT,NDG,IDGP,NDGP,NFLAG, &
                  NTRI,NCAVFC,IKEEP,NEDGRM,IEDGRM)
!
!     NUMBER NEW TETRAHEDRA REPLACING DELETED TETRAHEDRA WHERE
!     APPROPRIATE
!
      DO 40 K=1,NCNT
      IF (K.GT.NDEL) GO TO 35
      NEWCEL(K) = LDEL(K)
      GO TO 40
   35 NCELL     = NCELL  +1
      IF (NCELL.GT.MXCELL) GO TO 370
      NEWCEL(K) = NCELL
      NFLAG(NCELL) = 0
   40 CONTINUE
!
!     CHECK WHETHER A CAVITY BOUNDARY POINT IS BURIED IN THE CAVITY
!
      NTOT1     = 0
      DO 44 L=1,NDEL
      J         = LDEL(L)
      DO 44 I=1,4
      NCPNT     = NDC(I,J)
      IF (NFLAG(NCPNT).EQ.1) GO TO 44
      NTOT1     = NTOT1  +1
      IFLAG(NTOT1) = NCPNT
      NFLAG(NCPNT) = 1
   44 CONTINUE
      DO 46 K=1,NTOT1
      NFLAG(IFLAG(K)) = 0
   46 CONTINUE
      NTOT2     = 0
      DO 48 K=1,NCNT
      DO 48 I=1,3
      NCPNT     = NTRI(I,K)
      IF (NFLAG(NCPNT).EQ.1) GO TO 48
      NTOT2     = NTOT2  +1
      IFLAG(NTOT2) = NCPNT
      NFLAG(NCPNT) = 1
   48 CONTINUE
      DO 50 K=1,NTOT2
      NFLAG(IFLAG(K)) = 0
   50 CONTINUE
      IF (NTOT1.NE.NTOT2) GO TO 390
!
!     UPDATE NEIGHBORING TETRAHEDRON INFORMATION FOR UNDELETED
!     TETRAHEDRA
!
      DO 60 K=1,NCNT
      L         = NEWNBH(1,K)
      DO 55 M=1,4
      MM        = M
      IF (NBH(M,L).EQ.NOLD(K)) GO TO 60
   55 CONTINUE
      GO TO 380
   60 NSHAKE(K) = MM
!
!     GENERATE LIST OF CAVITY EDGES AND CONSTRUCT ARRAY NCAV
!
      CALL CAVEDG (NCNT,NEDG,IPOINT,NPOINT, &
                   NTRI,NCAV,NEWNBH,NEWCEL)
!
!     UPDATE FACE TO CELL POINTER FOR CAVITY FACES
!
      CALL DATSRF (NP,NCNT,NEDG,NDG,IDGP,NDGP,NEDGE, &
                   IPOINT,NCAV,NEDGRM,IEDGRM)
!
!     SET NBH ARRAY TO ZERO FOR THE DELETED TETRAHEDRA PRIOR TO
!     OVERWRITING WITH THE NEW TETRAHEDRON LIST
!
      DO 70 J=1,NDEL
      DO 70 I=1,4
   70 NBH(I,LDEL(J)) = 0
      DO 80 I=1,NEDG
      N1        = NCAV(1,I)
      N2        = NCAV(2,I)
      IPOINT(N1) = 0
      IPOINT(N2) = 0
   80 CONTINUE
!
!     COPY NEW TETRAHEDRON VERTEX AND NEIGHBOR INFORMATION INTO THE
!     ARRAYS NDC AND NBH
!
      DO 90 K=1,NCNT
      NBH(NSHAKE(K),NEWNBH(1,K)) = NEWCEL(K)
   90 CONTINUE
      DO 100 K=1,NCNT
      DO 95 I=1,4
      NBH(I,NEWCEL(K)) = NEWNBH(I,K)
   95 CONTINUE
      NDC(4,NEWCEL(K)) = NP
      NDC(1,NEWCEL(K)) = MIN(NTRI(1,K),NTRI(2,K),NTRI(3,K))
      NDC(2,NEWCEL(K)) = MAX(NTRI(1,K),NTRI(2,K),NTRI(3,K))
      NDC(3,NEWCEL(K)) = NTRI(1,K)  +NTRI(2,K)  +NTRI(3,K) &
                         -NDC(1,NEWCEL(K))  -NDC(2,NEWCEL(K))
      NPTET(NTRI(1,K)) = NEWCEL(K)
      NPTET(NTRI(2,K)) = NEWCEL(K)
      NPTET(NTRI(3,K)) = NEWCEL(K)
      XCEN(NEWCEL(K)) = XC(K)
      YCEN(NEWCEL(K)) = YC(K)
      ZCEN(NEWCEL(K)) = ZC(K)
      VOL(NEWCEL(K))  = V(K)
      RC(NEWCEL(K))   = RAD(K)
      RAT(NEWCEL(K))  = RCRIN(K)
      VOLMIN          = MIN(VOLMIN,V(K))
      RCMX            = MAX(RCMX,RAD(K))
  100 CONTINUE
      NPTET(NP) = NEWCEL(1)
!
!     COMPUTE CIRCUMRADIUS TO IN-RADIUS RATIO FOR NEW TETRAHEDRA
!
      DO 140 J=1,NDEL
      IF (NACPT(LDEL(J)).EQ.0) GO TO 140
      JPRE      = LNKUP(LDEL(J))
      JPAST     = LNKDN(LDEL(J))
      IF (JPRE.NE.0) LNKDN(JPRE) = JPAST
      IF (JPAST.NE.0) LNKUP(JPAST) = JPRE
      IF (JPRE.EQ.0) JFIRST = JPAST
      IF (JPAST.EQ.0) JLAST = JPRE
      NACPT(LDEL(J)) = 0
      NTRACK    = NTRACK  -1
  140 CONTINUE
      DO 150 K=1,NCNT
      N1       = NDC(1,NEWCEL(K))
      N2       = NDC(2,NEWCEL(K))
      N3       = NDC(3,NEWCEL(K))
      N4       = NDC(4,NEWCEL(K))
      NACPT(NEWCEL(K)) = 0
      DENT     = DENS(N1)  +DENS(N2)  +DENS(N3)  +DENS(N4)
      IF (RC(NEWCEL(K)).LT.0.22*DENT) GO TO 150
      NACPT(NEWCEL(K)) = 1
      NTRACK   = NTRACK  +1
      IF (JLAST.NE.0) LNKDN(JLAST) = NEWCEL(K)
      LNKDN(NEWCEL(K)) = 0
      LNKUP(NEWCEL(K)) = JLAST
      JLAST    = NEWCEL(K)
      IF (JFIRST.EQ.0) JFIRST = NEWCEL(K)
  150 CONTINUE
      RETURN
  295 CONTINUE
!     WRITE (6,595)
      DO 297 K=1,NDEL
      NFLAG(LDEL(K)) = 0
  297 CONTINUE
      NFAIL    = 5
      RETURN
  300 CONTINUE
      WRITE (6,600) NP
      WRITE (6,888) N1,X(1,N1),X(2,N1),X(3,N1), &
                    N2,X(1,N2),X(2,N2),X(3,N2), &
                    N3,X(1,N3),X(2,N3),X(3,N3), &
                    N4,X(1,N4),X(2,N4),X(3,N4),V(K),TOLV
  888 FORMAT('N1 ',I5,'   X = ',F8.4,'   Y = ',F8.4,'   Z = ',F8.4/ &
             'N2 ',I5,'   X = ',F8.4,'   Y = ',F8.4,'   Z = ',F8.4/ &
             'N3 ',I5,'   X = ',F8.4,'   Y = ',F8.4,'   Z = ',F8.4/ &
             'N4 ',I5,'   X = ',F8.4,'   Y = ',F8.4,'   Z = ',F8.4/ &
             'VOLUME = ',E13.5,' TOLV = ',E13.5)
      GO TO 312
  305 CONTINUE
!     WRITE (6,605) RADMAX,RAD(K),RATMAX,RCRIN(K)
      GO TO 312
  310 CONTINUE
      WRITE (6,610) V(K),TOLV
  312 DO 315 K=1,NDEL
      NFLAG(LDEL(K)) = 0
  315 CONTINUE
      NFAIL    = 1
      RETURN
  320 CONTINUE
!     WRITE (6,620)
      DO 325 K=1,NDEL
      NFLAG(LDEL(K)) = 0
  325 CONTINUE
      NFAIL    = 3
      RETURN
  350 CONTINUE
      VDIFF    = SUM1  -SUM2
      WRITE (6,650) VDIFF
      DO 355 K=1,NDEL
      NFLAG(LDEL(K)) = 0
  355 CONTINUE
      NFAIL    = 3
      RETURN
  360 CONTINUE
      WRITE (6,660)
      STOP
  370 CONTINUE
      WRITE (6,670)
      STOP
  380 CONTINUE
      WRITE (6,680)
      STOP
  390 CONTINUE
      WRITE (6,690) NTOT1,NTOT2
      STOP
  595 FORMAT(5X,'CAVITY DOES NOT CONTAIN ORIGINATING TETRAHEDRON')
  600 FORMAT(//5X,'AT LEAST ONE NEW TETRAHEDRON HAS TOO SMALL A VOLUME'/ &
               5X,'ADDRESS OF INSERTED POINT IS ',I6)
  605 FORMAT(/5X,'RADMAX ',E13.5,' RC ',E13.5, &
                 ' RATMAX ',E13.5,' RAT ',E13.5)
  610 FORMAT(5X,'VOLUME OF A NEW TETRAHEDRON IS LESS THAN TOLV'/ &
             5X,'VOLUME = ',E13.5,'  TOLV = ',E13.5/)
  620 FORMAT(5X,'NEW POINT CREATES A TETRAHEDRON TOO CLOSE TO BOUNDARY')
  650 FORMAT(5X,'VOLUME VISIBILITY CHECK FAILED, VDIFF = ',E13.5)
  670 FORMAT(//5X,'DIMENSION OF NDC EXCEEDED IN ROUTINE INSERT'/ &
               5X,'INCREASE SIZE OF MCELL.')
  660 FORMAT(//5X,'UNABLE TO FIND EDGE ADDRESS FOR A NEW TETRAHEDRON'/ &
               5X,'PROGRAM STOPPED IN ROUTINE INSERT')
  680 FORMAT(//5X,'UNABLE TO FIND A CONTIGUITY BETWEEN A NEW'/ &
               5X,'TETRAHEDRON AND A NON-CAVITY TETRAHEDRON'/ &
               5X,'PROGRAM STOPPED IN ROUTINE INSERT')
  690 FORMAT(//5X,'NUMBER OF CAVITY POINTS = ',I6,' IS DIFFERENT FROM'/ &
               5X,'THE NUMBER OF NEW CELL BOUNDARY POINTS = ',I6// &
               5X,'PROBABLE CAUSE IS A FAILURE IN THE TOLERANCE FOR'/ &
               5X,'THE DELAUNAY SPHERE TEST IN ROUTINE CAVITY.'/ &
               5X,'PROGRAM STOPPED IN ROUTINE INSERT.')
      END SUBROUTINE INSERT






!
!     ******************************************************************
!
      SUBROUTINE TETLOC (NP,NCLOSE,LBRK,LCONT,NFAIL, &
                         X,NDC,NBH,IPROT,DENS,VOL, &
                         IFLAG,NFLAG,NFILL,NEWCEL,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  STARTING WITH TETRAHEDRON LBRK, CARRY OUT A TREE SEARCH TO    *
!     *  TO FIND THE TETRAHEDRON LCONT THAT CONTAINS THE POINT NP.     *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      INTEGER :: LBRK,LCONT,NCLOSE,NFAIL,NP
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: NFILL(*),NEWCEL(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*),VOL(*)

      INTEGER :: K,K1,K2,L,LBEST,LCHK,LCNT,LNEXT,LSEEK,L1,L2,L3,L4, &
                 L5,NCHK,NCONT,NCNT,NMON,NPROX,N1,N2,N3,N4
      DOUBLE PRECISION :: DENB,RTEST,TOLS,TOLZ,VBEST,VDIFF,VTET, &
                          VTHRES,V1,V2,V3,V4,V5,XPT,YPT,ZPT
!
!     ******************************************************************
!
      TOLS      = 1.E-15
      TOLZ      = 1.E-8
      XPT       = X(1,NP)
      YPT       = X(2,NP)
      ZPT       = X(3,NP)
      LSEEK     = LBRK
      NMON      = 1
      IFLAG(NMON) = LSEEK

      CALL VOLCOM (LSEEK,NP,VDIFF,NCONT,X,NDC)

      VTHRES    = MAX(TOLV,TOLZ*VOL(LSEEK))
      IF (VDIFF.LT.VTHRES) GO TO 50
      IF (NCONT.EQ.0) GO TO 50
      IF (NCONT.EQ.-1) GO TO 360
      NPROX     = 0
      VBEST     = VDIFF
      LBEST     = LSEEK
      NFLAG(LSEEK) = 1
      NFILL(1)  = LSEEK
      LCNT      = 1
   10 NCNT      = 0
      NCHK      = 0
!
!     TREE SEARCH THROUGH NEIGHBORING TETRAHEDRA
!
      DO 20 LCHK=1,LCNT
      L         = NFILL(LCHK)
      DO 20 K=1,4
      LSEEK     = NBH(K,L)
      N1        = NDC(1,LSEEK)
      N2        = NDC(2,LSEEK)
      N3        = NDC(3,LSEEK)
      N4        = NDC(4,LSEEK)
      IF (NFLAG(LSEEK).EQ.1) GO TO 20
      IF (NDC(4,LSEEK).EQ.-1.AND.NPROX.EQ.0) GO TO 20
      IF (NDC(4,LSEEK).EQ.-1.AND.NPROX.GT.0) GO TO 15
      IF (NPROX.EQ.0.AND.N1.NE.NCLOSE.AND.N2.NE.NCLOSE &
                    .AND.N3.NE.NCLOSE.AND.N4.NE.NCLOSE) GO TO 20
      NCHK      = 1

      CALL VOLCOM (LSEEK,NP,VDIFF,NCONT,X,NDC)

      VTHRES    = MAX(TOLV,TOLZ*VOL(LSEEK))
      IF (VDIFF.LT.VTHRES) GO TO 50
      IF (NCONT.EQ.0) GO TO 50
      IF (NCONT.EQ.-1) GO TO 360
      IF (VDIFF.GT.VBEST) GO TO 15
      VBEST     = VDIFF
      LBEST     = LSEEK
   15 NCNT      = NCNT  +1
      IF (NCNT.GT.MXTEST) GO TO 310
      NEWCEL(NCNT) = LSEEK
      NMON      = NMON  +1
      IF (NMON.GT.MXNODE) GO TO 315
      IF (NMON.GT.MXNODE) GO TO 320
      IFLAG(NMON) = LSEEK
      NFLAG(LSEEK) = 1
   20 CONTINUE
      IF (NCNT.GT.0.AND.NPROX.EQ.0) GO TO 25
      IF (NPROX.EQ.0) GO TO 45
      IF (NPROX.EQ.15) GO TO 320
      NPROX     = NPROX  +1
      IF (NCNT.GT.0) GO TO 25

      RTEST     = VBEST/VOL(LBEST)
      IF (RTEST.GT.TOLZ) WRITE (6,920) NPROX,VBEST,VOL(LBEST),RTEST
  920 FORMAT('DIFFICULTY, NPROX = ',I4,' MIN VDIFF = ',E13.5, &
             ' VOL ',E13.5,' RTEST ',E13.5)

      IF (RTEST.GT.TOLZ) GO TO 320
      LSEEK     = LBEST
      GO TO 50
   25 LCNT      = NCNT
      DO 40 K=1,LCNT
      NFILL(K)  = NEWCEL(K)
   40 CONTINUE
      GO TO 10
   45 NPROX     = 1
      LSEEK     = LBRK
      NFILL(1)  = LSEEK
      LCNT      = 1 
      GO TO 10
!
!     TETRAHEDRON LSEEK CONTAINS POINT (XPT,YPT,ZPT).
!     INTERPOLATE LENGTH SCALE VALUE DENS FOR THE NEW POINT NP.
!
   50 DO 60 K=1,NMON
      NFLAG(IFLAG(K)) = 0
   60 CONTINUE
      LCONT     = LSEEK
      IF (IPROT(LCONT).EQ.1) GO TO 330
      L1        = NDC(1,LCONT)
      L2        = NDC(2,LCONT)
      L3        = NDC(3,LCONT)
      L4        = NDC(4,LCONT)

      CALL TETCOF (L1,L2,L3,L4,XPT,YPT,ZPT,V1,V2,V3,V4,VTET,X)

      DENS(NP)  = (V1*DENS(L1)  +V2*DENS(L2) &
                                +V3*DENS(L3)  +V4*DENS(L4))/VTET
      IF (V1.GE.TOLV.AND.V2.GE.TOLV.AND.V3.GE.TOLV &
                                   .AND.V4.GE.TOLV) GO TO 100
      IF (V2.LT.TOLS.AND.V3.LT.TOLS.AND.V4.LT.TOLS) GO TO 340
      IF (V1.LT.TOLS.AND.V3.LT.TOLS.AND.V4.LT.TOLS) GO TO 340
      IF (V1.LT.TOLS.AND.V2.LT.TOLS.AND.V4.LT.TOLS) GO TO 340
      IF (V1.LT.TOLS.AND.V2.LT.TOLS.AND.V3.LT.TOLS) GO TO 340
      IF (V1.LT.TOLV) GO TO 80
      IF (V2.LT.TOLV) GO TO 85
      IF (V3.LT.TOLV) GO TO 90

      CALL NEIGHB (LCONT,L1,L2,L3,LNEXT,K1,K2,NDC,NBH)

      L5        = NDC(1,LNEXT)  +NDC(2,LNEXT)  +NDC(3,LNEXT) &
                 +NDC(4,LNEXT)  -L1  -L2  -L3

      CALL TETCOF (L1,L2,L3,L5,XPT,YPT,ZPT,V1,V2,V3,V5,VTET,X)

      DENB      = (V1*DENS(L1)  +V2*DENS(L2)  +V3*DENS(L3) &
                                              +V5*DENS(L5))/VTET
      DENS(NP)  = .5*(DENS(NP)  +DENB)
      GO TO 100

   80 CALL NEIGHB (LCONT,L2,L3,L4,LNEXT,K1,K2,NDC,NBH)

      L5        = NDC(1,LNEXT)  +NDC(2,LNEXT)  +NDC(3,LNEXT) &
                 +NDC(4,LNEXT)  -L2  -L3  -L4

      CALL TETCOF (L2,L3,L4,L5,XPT,YPT,ZPT,V2,V3,V4,V5,VTET,X)

      DENB      = (V2*DENS(L2)  +V3*DENS(L3)  +V4*DENS(L4) &
                                              +V5*DENS(L5))/VTET
      DENS(NP)  = .5*(DENS(NP)  +DENB)
      GO TO 100

   85 CALL NEIGHB (LCONT,L3,L4,L1,LNEXT,K1,K2,NDC,NBH)

      L5        = NDC(1,LNEXT)  +NDC(2,LNEXT)  +NDC(3,LNEXT) &
                 +NDC(4,LNEXT)  -L3  -L4  -L1

      CALL TETCOF (L3,L4,L1,L5,XPT,YPT,ZPT,V3,V4,V1,V5,VTET,X)

      DENB      = (V3*DENS(L3)  +V4*DENS(L4)  +V1*DENS(L1) &
                                              +V5*DENS(L5))/VTET
      DENS(NP)  = .5*(DENS(NP)  +DENB)
      GO TO 100

   90 CALL NEIGHB (LCONT,L4,L1,L2,LNEXT,K1,K2,NDC,NBH)

      L5        = NDC(1,LNEXT)  +NDC(2,LNEXT)  +NDC(3,LNEXT) &
                 +NDC(4,LNEXT)  -L4  -L1  -L2

      CALL TETCOF (L4,L1,L2,L5,XPT,YPT,ZPT,V4,V1,V2,V5,VTET,X)

      DENB      = (V4*DENS(L4)  +V1*DENS(L1)  +V2*DENS(L2) &
                                              +V5*DENS(L5))/VTET
      DENS(NP)  = .5*(DENS(NP)  +DENB)
  100 RETURN
  310 CONTINUE
      WRITE (6,610)
      STOP
  315 CONTINUE
      WRITE (6,615)
      STOP
  320 CONTINUE
      DO 325 K=1,NMON
      NFLAG(IFLAG(K)) = 0
  325 CONTINUE
      NFAIL    = 1
      RETURN
  330 CONTINUE
!     WRITE (6,630)
      NFAIL    = 2
      RETURN
  340 CONTINUE
      WRITE (6,640) NP,X(1,NP),X(2,NP),X(3,NP), &
                    L1,L2,L3,L4,VTET,V1,V2,V3,V4
      STOP
  350 CONTINUE
      NFAIL    = 4
      RETURN
  360 CONTINUE
      WRITE (6,660)
 
      WRITE (6,910) NP,X(1,NP),X(2,NP),X(3,NP),NCONT
  910 FORMAT('NP = ',I6,' X = ',F8.3,' Y = ',F8.3,' Z = ',F8.3, &
             ' NCONT ',I2)
      STOP
 
  610 FORMAT(//5X,'DIMENSION OF NEWCEL EXCEEDED IN ROUTINE TETLOC'/ &
               5X,'INCREASE SIZE OF MTEST')
  615 FORMAT(//5X,'DIMENSION OF IFLAG EXCEEDED IN ROUTINE TETLOC'/ &
               5X,'INCREASE SIZE OF MBPTS')
  620 FORMAT(5X,'UNABLE TO FIND TETRAHEDRON THAT CONTAINS NEW POINT', &
                ' IN ROUTINE TETLOC')
  630 FORMAT(5X,'POINT LIES OUTSIDE DOMAIN TO BE MESHED')
  640 FORMAT(//5X,'NEW POINT, NP = ',I6,' X = ',F6.2,' Y = ',F6.2, &
             ' Z = ',F6.2/5X,'APPEARS TO BE COINCIDENT WITH AN', &
             ' EXISTING POINT IN THE MESH'/ &
             5X,'VERTEX ADDRESSES OF CONTAINING TETRAHEDRON'/ &
             5X,'L1 ',I6,' L2 ',I6,' L3 ',I6,' L4 ',I6,' VOL = ',E13.5/ &
             5X,'V1 = ',E13.5,' V2 = ',E13.5,' V3 = ',E13.5, &
             ' V4 = ',E13.5/5X,'PROGRAM STOPPED IN ROUTINE TETLOC')
  660 FORMAT(//5X,'A GEOMETRIC INCONSISTENCY HAS BEEN FOUND IN ROUTINE', &
                  '  TETLOC.'/5X,'THE INSERTED POINT HAS A NEGATIVE', &
                  ' ORIENTATION WITH RESPECT TO EVERY TETRAHEDRON FACE')
      END SUBROUTINE TETLOC






!
!     ******************************************************************
!
      SUBROUTINE CAVITY (NP,LCONT,NDEL,LDEL, &
                         X,NDC,NBH,IPROT,IFLAG,NFLAG, &
                         XCEN,YCEN,ZCEN,RC,NFILL,NEWCEL,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  STARTING FROM TETRAHEDRON LCONT WHICH CONTAINS THE POINT NP,  *
!     *  CARRY OUT A TREE SEARCH TO FIND THE COMPLETE CAVITY OF        *
!     *  TETRAHEDRA WHOSE CIRCUMSPHERES CONTAIN POINT NP               *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: LCONT,NDEL,NP
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: NFILL(*),NEWCEL(*)
      INTEGER :: LDEL(*)
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: XCEN(*),YCEN(*),ZCEN(*),RC(*)
      
      INTEGER :: L,LSEEK,L1,L2,L3,K,K1,K2,LCHK,LCNT,LNEXT,NCNT,NMON,M, &
                 N1,N2,N3,N4,N5
      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,DCEN,RNX,RNY,RNZ, &
                          THRES,TOLP,VTET1,VTET2,XPT,YPT,ZPT     
!
!     ******************************************************************
!
      TOLP      = 1.E-12
      XPT       = X(1,NP)
      YPT       = X(2,NP)
      ZPT       = X(3,NP)
      NDEL      = 1
      LDEL(1)   = LCONT
      LCNT      = 1
      NFILL(1)  = LCONT
      NMON      = 1
      IFLAG(1)  = LCONT
      NFLAG(LCONT) = 1
   10 NCNT      = 0
!
!     SEARCH THROUGH NEIGHBORING TETRAHEDRA
!
      DO 20 LCHK=1,LCNT
      L         = NFILL(LCHK)
      DO 20 K=1,4
      LSEEK     = NBH(K,L)
!
!     CHECK WHETHER NEIGHBORING TETRAHEDRON LSEEK HAS ALREADY
!     BEEN EXAMINED
!
      IF (NFLAG(LSEEK).EQ.1) GO TO 20
!
!     CHECK WHETHER NEIGHBORING TETRAHEDRON LSEEK LIES OUTSIDE
!     THE CONVEX HULL
!
      IF (NDC(4,LSEEK).EQ.-1) GO TO 20
      NMON      = NMON  +1
      IF (NMON.GT.MXNODE) GO TO 315
      IFLAG(NMON) = LSEEK
      NFLAG(LSEEK) = 1
!
!     CHECK WHETHER CIRCUMSPHERE OF TETRAHEDRON LSEEK CONTAINS
!     THE POINT (XPT,YPT,ZPT)
!
      DCEN      = SQRT((XPT  -XCEN(LSEEK))**2  +(YPT  -YCEN(LSEEK))**2 &
                                               +(ZPT  -ZCEN(LSEEK))**2)
      IF (DCEN.GE.RC(LSEEK)*(1.+TOLP))GO TO 20
!
!     CHECK WHETHER COMMON FACE IS PROTECTED
!
      IF (IPROT(LSEEK).EQ.1) GO TO 20
!
!     CHECK WHETHER ANY FACES OF TETRAHEDRON LSEEK ARE PROTECTED.
!     IF A FACE IS PROTECTED, CHECK WHETHER FACE IS VISIBLE
!     FROM POINT (XPT,YPT,ZPT)
!
      DO 15 M=1,4
      LNEXT     = NBH(M,LSEEK)
      IF (IPROT(LNEXT).EQ.0) GO TO 15

      CALL LOCK (LSEEK,LNEXT,N1,N2,N3,N4,N5,K1,K2,NDC,NBH)

      L1      = MAX(N1,N2,N3)
      L3      = MIN(N1,N2,N3)
      L2      = N1  +N2  +N3  -L1  -L3
      RNX     = COFACT (X(2,L1),X(2,L2),X(2,L3),X(3,L1),X(3,L2),X(3,L3))
      RNY     = COFACT (X(3,L1),X(3,L2),X(3,L3),X(1,L1),X(1,L2),X(1,L3))
      RNZ     = COFACT (X(1,L1),X(1,L2),X(1,L3),X(2,L1),X(2,L2),X(2,L3))
      VTET1   = RNX*(XPT  -X(1,L1))  +RNY*(YPT  -X(2,L1)) &
                                     +RNZ*(ZPT  -X(3,L1))
      VTET2   = RNX*(X(1,N4)  -X(1,L1))  +RNY*(X(2,N4)  -X(2,L1)) &
                                         +RNZ*(X(3,N4)  -X(3,L1))
      THRES   = MAX(1.0D0,ABS(VTET2))
      IF (ABS(VTET1).LT.TOLV*THRES) GO TO 20
      IF (VTET1*VTET2.LT.0.) GO TO 20
   15 CONTINUE
!
!     ADMIT TETRAHEDRON LSEEK TO LIST OF CAVITY TETRAHEDRA
!
   17 NCNT      = NCNT  +1
      IF (NCNT.GT.MXTEST) GO TO 310
      NEWCEL(NCNT) = LSEEK
      NDEL      = NDEL  +1
      IF (NDEL.GT.MXTEST) GO TO 320
      LDEL(NDEL) = LSEEK
   20 CONTINUE
      IF (NCNT.EQ.0) GO TO 40
      LCNT      = NCNT
      DO 30 K=1,LCNT
      NFILL(K)  = NEWCEL(K)
   30 CONTINUE
      GO TO 10
   40 DO 50 K=1,NMON
      NFLAG(IFLAG(K)) = 0
   50 CONTINUE
      DO 60 K=1,NDEL
      NFLAG(LDEL(K)) = 1
   60 CONTINUE
      RETURN
  310 CONTINUE
      WRITE (6,610)
      STOP
  315 CONTINUE
      WRITE (6,615)
      STOP
  320 CONTINUE
      WRITE (6,620)
      STOP
  610 FORMAT(//5X,'DIMENSION OF NEWCEL EXCEEDED IN ROUTINE CAVITY'/ &
               5X,'INCREASE SIZE OF MTEST')
  615 FORMAT(//5X,'DIMENSION OF IFLAG EXCEEDED IN ROUTINE CAVITY'/ &
               5X,'INCREASE SIZE OF MNODE')
  620 FORMAT(//5X,'DIMENSION OF LDEL EXCEEDED IN ROUTINE CAVITY'/ &
               5X,'INCREASE SIZE OF MTEST')
      END SUBROUTINE CAVITY




!
!     ******************************************************************
!
      SUBROUTINE CAVBND (NP,LCONT,NDEL,LDEL,NCNT, &
                         X,NDC,NBH,VOL,IFLAG,NFLAG, &
                         NFILL,NEWCEL,NTRI,NEWNBH,NOLD,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  SEARCH THROUGH LIST OF CAVITY TETRAHEDRA TO FIND THE FACES    *
!     *  NTRI(3,--) ON THE CAVITY BOUNDARY                             *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: LCONT,NCNT,NDEL,NP
      INTEGER :: NDC(4,*),NBH(4,*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: NTRI(3,*),NFILL(*),NEWNBH(4,*),NOLD(*),NEWCEL(*)
      INTEGER :: LDEL(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*),VOL(*)
          
      INTEGER :: J,K,L,LCHK,LCNT,LL,LTRY,M,MEND,MMAX,MMIN,MSUM, &
                 M1,M2,M3,M4,N,NCONT,NEND,NJOIN,NMAX,NMON,NMIN, &
                 NSUM,NWATCH,N1,N2,N3,N4
      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,FAC,RNX,RNY,RNZ,VDIFF, &
                          VTET,VTET1,VTET2,V1,V2,V3,V4,XPT,YPT,ZPT      
!
!     ******************************************************************
!
      XPT       = X(1,NP)
      YPT       = X(2,NP)
      ZPT       = X(3,NP)
   10 NCNT      = 0
      NWATCH    = 0
      DO 70 K=1,NDEL
      NJOIN     = 0
      N1        = NDC(1,LDEL(K))
      N2        = NDC(2,LDEL(K))
      N3        = NDC(3,LDEL(K))
      N4        = NDC(4,LDEL(K))
      DO 65 J=1,4
      L         = NBH(J,LDEL(K))
      IF (NFLAG(L).GT.0) GO TO 65
      NEND      = N4
      M1        = NDC(1,L)
      M2        = NDC(2,L)
      M3        = NDC(3,L)
      M4        = NDC(4,L)
   20 MEND      = M4
      NMIN      = MIN(N1,N2,N3)
      NMAX      = MAX(N1,N2,N3)
      NSUM      = N1  +N2  +N3
   30 MMIN      = MIN(M1,M2,M3)
      MMAX      = MAX(M1,M2,M3)
      MSUM      = M1  +M2  +M3
      IF (MMIN.EQ.NMIN.AND.MMAX.EQ.NMAX.AND.MSUM.EQ.NSUM) GO TO 40
      IF (M1.EQ.MEND) GO TO 35
      M         = M1
      M1        = M2
      M2        = M3
      M3        = M4
      M4        = M
      GO TO 30
   35 IF (N1.EQ.NEND) GO TO 300
      N         = N1
      N1        = N2
      N2        = N3
      N3        = N4
      N4        = N
      GO TO 20
!
!     CHECK WHETHER CAVITY FACE (N1,N2,N3) IS VISIBLE FROM 
!     POINT (XPT,YPT,ZPT)
!
   40 RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      FAC     = 1./SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      VTET1   = FAC*(RNX*(XPT  -X(1,N1))  +RNY*(YPT  -X(2,N1)) &
                                          +RNZ*(ZPT  -X(3,N1)))
      VTET2   = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                         +RNZ*(X(3,N4)  -X(3,N1))
!     THRES   = MAX(1.,ABS(VTET2))
!     IF (LDEL(K).NE.LCONT.AND.ABS(VTET1).LT.TOLV*THRES) GO TO 45
      IF (LDEL(K).NE.LCONT.AND.ABS(VTET1).LT.TOLV) GO TO 45
      IF (VTET1*VTET2.GE.0.) GO TO 50
   45 NFLAG(LDEL(K)) = 2
      NWATCH  = 1
      GO TO 70
!
!     CAVITY FACE (N1,N2,N3) IS VISIBLE. STORE VERTICES IN ARRAY
!     NTRI
!
   50 NCNT      = NCNT  +1
      IF (NCNT.GT.MXTEST) GO TO 310
      NJOIN     = NJOIN  +1
      NTRI(1,NCNT) = N1
      NTRI(2,NCNT) = N2
      NTRI(3,NCNT) = N3
      NEWNBH(1,NCNT) = L
      NOLD(NCNT)  = LDEL(K)
   65 CONTINUE
      IF (NJOIN.EQ.4.AND.NDEL.GT.1) GO TO 320
   70 CONTINUE
      IF (NWATCH.EQ.0) RETURN
      IF (NFLAG(LCONT).EQ.2) GO TO 330
!
!     REPEAT TREE SEARCH TO FIND A CONTIGUOUS SET OF TETRAHEDRA,
!     STARTING FROM TETRAHEDRON THAT CONTAINS POINT NP
!
      LCNT      = 1
      NFILL(1)  = LCONT
      NMON      = 1
      IFLAG(1)  = LCONT
      NFLAG(LCONT) = 2
   80 NCNT      = 0
      DO 85 LCHK=1,LCNT
      L         = NFILL(LCHK)
      DO 85 K=1,4
      LTRY      = NBH(K,L)
      IF (NFLAG(LTRY).NE.1) GO TO 85
      NFLAG(LTRY) = 2
      NCNT      = NCNT  +1
      NEWCEL(NCNT) = LTRY
      NMON      = NMON  +1
      IF (NMON.GT.MXNODE) GO TO 335
      IFLAG(NMON) = LTRY
   85 CONTINUE
      IF (NCNT.EQ.0) GO TO 95
      LCNT      = NCNT
      DO 90 K=1,NCNT
      NFILL(K)  = NEWCEL(K)
   90 CONTINUE
      GO TO 80
!
!     RE-INITIALIZE NFLAG ARRAY
!
   95 DO 100 K=1,NDEL
      NFLAG(LDEL(K)) = 0
  100 CONTINUE
      NDEL      = NMON
      DO 110 K=1,NDEL
      LDEL(K)   = IFLAG(K)
      NFLAG(LDEL(K)) = 1
  110 CONTINUE
      GO TO 10
  300 CONTINUE
      WRITE (6,600)
      STOP
  310 CONTINUE
      WRITE (6,610)
      STOP
  320 CONTINUE
      WRITE (6,620)
      STOP
  330 CONTINUE
      WRITE (6,630)

      CALL VOLCOM (LCONT,NP,VDIFF,NCONT,X,NDC)

      N1      = NDC(1,LCONT)
      N2      = NDC(2,LCONT)
      N3      = NDC(3,LCONT)
      N4      = NDC(4,LCONT)

      CALL TETCOF (N1,N2,N3,N4,XPT,YPT,ZPT,V1,V2,V3,V4,VTET,X)

      WRITE (6,998) NP,NDEL,LCONT,N1,N2,N3,N4, &
                    V1,V2,V3,V4,VTET,VDIFF,NCONT
  998 FORMAT('NP ',I6,' NDEL ',I4,' LDEL ',I7,' N1,N2,N3,N4 ',4I6/ &
       'V1 ',E13.5,' V2 ',E13.5,' V3 ',E13.5,' V4 ',E13.5/ &
       ' VTET ',E13.5,' VDIFF ',E13.5,' NCONT ',I2)
      DO 996 K=1,NDEL
      LL      = LDEL(K)
      N1      = NDC(1,LL)
      N2      = NDC(2,LL)
      N3      = NDC(3,LL)
      N4      = NDC(4,LL)
  996 WRITE (6,997) LL,NFLAG(LL),N1,N2,N3,N4,VOL(LL)
  997 FORMAT('CELL ',I7,' NFLAG ',I2,' N1,N2,N3,N4 ',4I6,' VOL ',E13.5)
      STOP
  335 CONTINUE
      WRITE (6,635)
      STOP
  600 FORMAT(//5X,'UNABLE TO FIND A CONTIGUITY BETWEEN A CAVITY'/ &
               5X,'TETRAHEDRON AND A NON-CAVITY TETRAHEDRON'/ &
               5X,'PROGRAM STOPPED IN ROUTINE CAVBND')
  610 FORMAT(//5X,'DIMENSION OF NTRI ARRAY EXCEEDED IN ROUTINE CAVBND'/ &
               5X,'INCREASE SIZE OF MTEST')
  620 FORMAT(//5X,'FOUR CONTIGUITIES HAVE BEEN FOUND FOR A NEW' &
               ,' TETRAHEDRON'/ &
               5X,'NEW TETRAHEDRON IS THEREFORE ISOLATED'/ &
               5X,'PROGRAM STOPPED IN ROUTINE CAVBND') 
  630 FORMAT(//5X,'TETRAHEDRON CONTAINING POINT HAS FAILED'/ &
               5X,'VISIBILITY TEST'/ &
               5X,'PROGRAM STOPPED IN ROUTINE CAVBND')
  635 FORMAT(//5X,'DIMENSION OF IFLAG EXCEEDED IN ROUTINE CAVBND.'/ &
               5X,'INCREASE SIZE OF MBPTS')
      END SUBROUTINE CAVBND






!
!     ******************************************************************
!
      SUBROUTINE RECON (LDEL,NDEL,NCNT, &
                        NDC,NBH,IPROT,NDG,IDGP,NDGP,NFLAG, &
                        NTRI,NCAVFC,IKEEP,NEDGRM,IEDGRM)
!
!     ******************************************************************
!     *                                                                *
!     *  REMOVE INTERNAL CAVITY EDGES FROM DATA STRUCTURE NDG.         *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IEDGRM,NCNT,NDEL
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*),NFLAG(*)
      INTEGER :: NTRI(3,*),NEDGRM(*)
      INTEGER :: IDGP(*),NDGP(*)
      INTEGER :: LDEL(*),NCAVFC(3,*),IKEEP(*)

      INTEGER :: J,K,KCNT,L,L1,L2,M,MEND,MMAX,MMIN,MSUM,M1,M2,M3,M4, &
                 N,NEDG,NEND,NF,NMAX,NMIN,NSUM,N1,N2,N3,N4
!
!     ******************************************************************
!
!     FORM LIST OF INTERNAL CAVITY FACES
!
      IEDGRM    = 0
      NF        = 0
      DO 70 K=1,NDEL
      IF (IPROT(LDEL(K)).EQ.1) GO TO 370
      NFLAG(LDEL(K)) = 0
      DO 65 J=1,4
      L         = NBH(J,LDEL(K))
      IF (NFLAG(L).EQ.0) GO TO 65
      N1        = NDC(1,LDEL(K))
      N2        = NDC(2,LDEL(K))
      N3        = NDC(3,LDEL(K))
      N4        = NDC(4,LDEL(K))
      NEND      = N4
      M1        = NDC(1,L)
      M2        = NDC(2,L)
      M3        = NDC(3,L)
      M4        = NDC(4,L)
   40 MEND      = M4
      NMIN      = MIN(N1,N2,N3)
      NMAX      = MAX(N1,N2,N3)
      NSUM      = N1  +N2  +N3
   45 MMIN      = MIN(M1,M2,M3)
      MMAX      = MAX(M1,M2,M3)
      MSUM      = M1  +M2  +M3
      IF (MMIN.EQ.NMIN.AND.MMAX.EQ.NMAX.AND.MSUM.EQ.NSUM) GO TO 60
      IF (M1.EQ.MEND) GO TO 50
      M         = M1
      M1        = M2
      M2        = M3
      M3        = M4
      M4        = M
      GO TO 45
   50 IF (N1.EQ.NEND) GO TO 300
      N         = N1
      N1        = N2
      N2        = N3
      N3        = N4
      N4        = N
      GO TO 40
   60 NF        = NF  +1
      IF (NF.GT.MXTEST) GO TO 310
      NCAVFC(1,NF) = N1
      NCAVFC(2,NF) = N2
      NCAVFC(3,NF) = N3
   65 CONTINUE
   70 CONTINUE
      IF (NF.EQ.0) RETURN
!
!     FORM LIST OF EDGES ON BOUNDARY CAVITY
!
      KCNT      = 0
      DO 100 L=1,NCNT
      N1        = NTRI(1,L)
      N2        = NTRI(2,L)
      N3        = NTRI(3,L)
   90 L1        = MIN(N1,N2)
      L2        = MAX(N1,N2)
      NEDG       = IDGP(L1)
   95 IF (NEDG.EQ.0) GO TO 320
      IF (L2.EQ.NDG(2,NEDG)) GO TO 97
      NEDG       = NDGP(NEDG)
      GO TO  95
   97 IF (NFLAG(NEDG).EQ.1) GO TO 98
      NFLAG(NEDG) = 1
      KCNT      = KCNT  +1
      IF (KCNT.GT.MXTEST) GO TO 330
      IKEEP(KCNT) = NEDG
   98 IF (N1.EQ.NTRI(3,L)) GO TO 100
      N         = N1
      N1        = N2
      N2        = N3
      N3        = N
      GO TO 90
  100 CONTINUE
!
!     CHECK FOR CAVITY EDGES THAT ARE NOT ON BOUNDARY
!
      DO 140 J=1,NF
      N1        = NCAVFC(1,J)
      N2        = NCAVFC(2,J)
      N3        = NCAVFC(3,J)
  110 L1        = MIN(N1,N2)
      L2        = MAX(N1,N2)
      NEDG       = IDGP(L1)
  120 IF (NEDG.EQ.0) GO TO 320
      IF (L2.EQ.NDG(2,NEDG)) GO TO 130
      NEDG       = NDGP(NEDG)
      GO TO 120
  130 IF (NFLAG(NEDG).NE.0) GO TO 135
      NFLAG(NEDG) = 2
      KCNT      = KCNT  +1
      IF (KCNT.GT.MXTEST) GO TO 330
      IKEEP(KCNT) = NEDG
  135 IF (N1.EQ.NCAVFC(3,J)) GO TO 140
      N         = N1
      N1        = N2
      N2        = N3
      N3        = N
      GO TO 110
  140 CONTINUE
      IF (KCNT.EQ.0) RETURN
      DO 180 K=1,KCNT
      NEDG      = IKEEP(K)
      IF (NFLAG(NEDG).NE.2) GO TO 170

      CALL EDGERM (NEDG,NDG,IDGP,NDGP)

      IF (IEDGRM.EQ.MXTEST) GO TO 170
      IEDGRM    = IEDGRM  +1
      NEDGRM(IEDGRM) = NEDG
  170 NFLAG(NEDG) = 0
  180 CONTINUE
      RETURN

  300 WRITE (6,600)
      STOP
  310 WRITE (6,610)
      STOP
  320 WRITE (6,620)
      STOP
  330 WRITE (6,630)
      STOP
  350 WRITE (6,650)
      STOP
  370 WRITE (6,670)
      STOP
  380 WRITE (6,680)
      STOP
  600 FORMAT(//5X,'UNABLE TO FIND A CONTIGUITY BETWEEN CAVITY'/ &
               5X,'TETRAHEDRA'/ &
               5X,'PROGRAM STOPPED IN ROUTINE RECON') 
  610 FORMAT(//5X,'DIMENSION OF NTRI ARRAY EXCEEDED IN ROUTINE RECON'/ &
               5X,'INCREASE SIZE OF MTEST')
  620 FORMAT(5X,'UNABLE TO FIND EDGE IN CAVITY AMONG NDG ARRAY',/ &
                'PROGRAM STOPPED IN ROUTINE RECON')
  630 FORMAT(//5X,'DIMENSION OF IKEEP ARRAY EXCEEDED IN ROUTINE RECON'/ &
               5X,'INCREASE SIZE OF MTEST')
  650 FORMAT(5X,'AT LEAST ONE INTERNAL CAVITY FACE IS A PROTECTED'/ &
             5X,'BOUNDARY FACE. PROGRAM STOPPED IN ROUTINE RECON')
  670 FORMAT(/5X,'AT LEAST ONE CAVITY TETRAHEDRON IS PROTECTED.'/ &
              5X,'PROGRAM STOPPED IN ROUTINE RECON.')
  680 FORMAT(//5X,'UNABLE TO FIND FACE IN LINKED LIST'/ &
               5X,'PROGRAM STOPPED IN RECON')
      END SUBROUTINE RECON






!
!     ******************************************************************
!
      SUBROUTINE CAVEDG (NCNT,NEDG,IPOINT,NPOINT, &
                         NTRI,NCAV,NEWNBH,NEWCEL)
!
!     ******************************************************************
!     *                                                                *
!     *  GENERATE LIST OF EDGES ON CAVITY BOUNDARY AND CONSTRUCT EDGE  *
!     *  ARRAY NCAV(4,--)                                              *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: NCNT,NEDG
      INTEGER :: NPOINT(*),IPOINT(*)
      INTEGER :: NTRI(3,*),NEWNBH(4,*),NEWCEL(*),NCAV(4,*)

      INTEGER :: I,J1,J2,K,L,L1,L2,L3,NEXT,N1,N2
!
!     ******************************************************************
!
      NEDG      = 0
      DO 60 K=1,NCNT
      L1        = NTRI(1,K)
      L2        = NTRI(2,K)
      L3        = NTRI(3,K)
   10 N1        = MAX(L1,L2)
      N2        = MIN(L1,L2)
      I         = IPOINT(N1)
      IF (I.EQ.0) GO TO 30
   20 IF (N2.EQ.NCAV(2,I)) GO TO 40
      NEXT      = NPOINT(I)
      IF (NEXT.EQ.0) GO TO 30
      I         = NEXT
      GO TO 20
   30 NEDG      = NEDG  +1
      IF (NEDG.GT.MXCAV) GO TO 300
      NCAV(1,NEDG) = N1
      NCAV(2,NEDG) = N2
      NCAV(3,NEDG) = K
      NCAV(4,NEDG) = 0
      NPOINT(NEDG) = 0
      IF (I.NE.0) NPOINT(I)  = NEDG
      IF (I.EQ.0) IPOINT(N1) = NEDG
      GO TO 50
   40 NCAV(4,I) = K
   50 IF (L1.EQ.NTRI(3,K)) GO TO 60
      L         = L1
      L1        = L2
      L2        = L3
      L3        = L
      GO TO 10
   60 CONTINUE
!
!     FROM EDGE LIST,FIND REMAINING THREE NEIGHBORING
!     TETRAHEDRA TO EACH NEW VERTEX
!
      DO 70 K=1,NCNT
   70 NPOINT(K) = 1
      DO 80 I=1,NEDG
      J1        = NCAV(3,I)
      J2        = NCAV(4,I)
      NPOINT(J1) = NPOINT(J1)  +1
      NPOINT(J2) = NPOINT(J2)  +1
      NEWNBH(NPOINT(J1),J1) = NEWCEL(J2)
      NEWNBH(NPOINT(J2),J2) = NEWCEL(J1)
   80 CONTINUE
      DO 90 K=1,NCNT
      IF (NPOINT(K).NE.4) GO TO 310
   90 CONTINUE
!
!     SET IPOINT ENTRIES BACK TO ZERO
!
      DO 100 K=1,NCNT
      IF (NTRI(1,K).GT.0) IPOINT(NTRI(1,K)) = 0
      IF (NTRI(2,K).GT.0) IPOINT(NTRI(2,K)) = 0
      IF (NTRI(3,K).GT.0) IPOINT(NTRI(3,K)) = 0
  100 CONTINUE
      RETURN
  300 CONTINUE
      WRITE (6,600)
      STOP
  310 CONTINUE
      WRITE (6,610)
      STOP
  600 FORMAT(//5X,'DIMENSION OF NCAV EXCEEDED IN ROUTINE CAVEDG'/ &
               5X,'INCREASE SIZE OF MCAV')
  610 FORMAT(//5X,'UNABLE TO FIND ALL CONTIGUITIES BETWEEN NEW', &
                  ' TETRAHEDRA'/5X,'PROGRAM STOPPED IN ROUTINE CAVEDG')
      END SUBROUTINE CAVEDG





!
!     ******************************************************************
!
      SUBROUTINE DATSRF (NP,NCNT,NEDG,NDG,IDGP,NDGP,NEDGE, &
                         IPOINT,NCAV,NEDGRM,IEDGRM)
!
!     ******************************************************************
!     *                                                                *
!     *  UPDATE EDGE DATA STRUCTURE NDG                                *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IEDGRM,NCNT,NEDG,NEDGE,NP
      INTEGER :: NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: IPOINT(*)
      INTEGER :: NCAV(4,*),NEDGRM(*)

      INTEGER :: I,IEDG,J1,J2,L1,L2,LP,NEDCNT,NEWED,N1,N2

!
!     ******************************************************************
!
!     INSERT NEW EDGES INTO LIST NDG
!
      NEDCNT    = 0
      DO 80 I=1,NEDG
      N1        = NCAV(1,I)
      N2        = NCAV(2,I)
      J1        = NCAV(3,I)
      J2        = NCAV(4,I)
      IF (IPOINT(N1).GT.0) GO TO 65
      L1        = MIN(N1,NP)
      LP        = MAX(N1,NP)
      IEDG      = IDGP(L1)
      IF (IEDG.EQ.0) GO TO 50
   40 IF (NDGP(IEDG).EQ.0) GO TO 50
      IEDG      = NDGP(IEDG)
      GO TO 40
   50 NEDCNT    = NEDCNT  +1
      IF (NEDCNT.GT.IEDGRM) GO TO 55
      NEWED     = NEDGRM(NEDCNT)
      GO TO 60
   55 NEDGE     = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 310
      NEWED     = NEDGE
   60 NDG(1,NEWED) = L1
      NDG(2,NEWED) = LP
      NDGP(NEWED) = 0
      IF (IEDG.EQ.0) IDGP(L1) = NEWED
      IF (IEDG.GT.0) NDGP(IEDG) = NEWED
      IPOINT(N1) = 1
   65 IF (IPOINT(N2).GT.0) GO TO 80
      L2        = MIN(N2,NP)
      LP        = MAX(N2,NP)
      IEDG      = IDGP(L2)
      IF (IEDG.EQ.0) GO TO 68
   66 IF (NDGP(IEDG).EQ.0) GO TO 68
      IEDG      = NDGP(IEDG)
      GO TO 66
   68 NEDCNT    = NEDCNT  +1
      IF (NEDCNT.GT.IEDGRM) GO TO 70
      NEWED     = NEDGRM(NEDCNT)
      GO TO 75
   70 NEDGE     = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 310
      NEWED     = NEDGE
   75 NDG(1,NEWED) = L2
      NDG(2,NEWED) = LP
      NDGP(NEWED) = 0
      IF (IEDG.EQ.0) IDGP(L2) = NEWED
      IF (IEDG.GT.0) NDGP(IEDG) = NEWED
      IPOINT(N2) = 1
   80 CONTINUE
      RETURN
  310 WRITE (6,610)
      STOP
  610 FORMAT(//5X,'DIMENSION OF NDG ARRAY EXCEEDED IN ROUTINE DATSRF.'/ &
               5X,'INCREASE SIZE OF MEDGE.')
      END SUBROUTINE DATSRF





!
!     ******************************************************************
!
      SUBROUTINE COLAPS (NA,NB,NCOL,NVERT,NFAIL, &
                         X,NDC,NBH,IPROT, &
                         ITYP,XCEN,YCEN,ZCEN,VOL,RC,RAT, &
                         NVCNT,IFLAG,NFLAG,NPTET, &
                         NDG,IDGP,NDGP, &
                         NOCTR,IOCTR,NLINK,XFAR,YFAR,ZFAR,IDONE,NREF, &
                         KSRCH,NSRCH,IRING,NTETKP, &
                         LNBR,ISHK,MNBR,KSHK,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  COLLAPSE EDGE (NA,NB) AND REMOVE THE TETRAHEDRA SURROUNDING   *
!     *  EACH COLLAPSED EDGE.                                          *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1998                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,NA,NB,NCOL,NFAIL,NVERT
      INTEGER :: ITYP(*),NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*)
      INTEGER :: IDGP(*),NDGP(*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: NVCNT(*),NPTET(*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: NSRCH(*),KSRCH(*)
      INTEGER :: IRING(*),NTETKP(*),LNBR(*),ISHK(*),MNBR(*),KSHK(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)

      INTEGER :: I,IEDG,IMON,INAD,IP1,ISMALL,ITYPT,JNEXT,JPRE,K,KCNT, &
                 KK,K1,K2,L,LC,LN,L1,L2,M,M1,M2,M3,M4,N,NABDY,NBBDY, &
                 NCNT,NEXCH,NKEEP,NP,NPASS,NPNEXT,NPRE,NRBDY,NRNGP1, &
                 NRING,N1,N2,N3,N4
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,RADC,RADMAX,RATMAX, &
                          RCRIN,RNX,RNY,RNZ,TOLS,VL,VTET1,VTET2,XCN, &
                          YCN,ZCN,XPT,YPT,ZPT
!
!     ******************************************************************
!
      NFAIL     = 0
      NVERT     = 0
      TOLS      = 1.E-9
      IF (NA.LT.0.OR.NB.LT.0) RETURN
      IF (ITYP(NA).LT.0.OR.ITYP(NB).LT.0) RETURN
!
!     FIND A TETRAHEDRON THAT IS INCIDENT TO EDGE (NA,NB)
!
      KSRCH(1)  = NPTET(NA)
      IMON      = 0
      KCNT      = 1
   10 NCNT      = 0
      DO 22 K=1,KCNT
      LC        = KSRCH(K)
      IF (LC.LE.0) GO TO 300
      IF (NFLAG(LC).EQ.1) GO TO 22
      M1        = NDC(1,LC)
      M2        = NDC(2,LC)
      M3        = NDC(3,LC)
      M4        = NDC(4,LC)
      IF (M1.NE.NB.AND.M2.NE.NB.AND.M3.NE.NB.AND.M4.NE.NB) GO TO 15
      GO TO 30
   15 NFLAG(LC)  = 1
      IMON      = IMON  +1
      IF (IMON.GT.MXNODE) GO TO 310
!     IF (IMON.GT.MXNODE) GO TO 200
      IFLAG(IMON) = LC
      DO 20 M=1,4
      LN        = NBH(M,LC)
      IF (NFLAG(LN).EQ.1) GO TO 20
      M1        = NDC(1,LN)
      M2        = NDC(2,LN)
      M3        = NDC(3,LN)
      M4        = NDC(4,LN)
      IF (M1.NE.NA.AND.M2.NE.NA.AND.M3.NE.NA.AND.M4.NE.NA) GO TO 20
      NCNT      = NCNT  +1
      IF (NCNT.GT.MXTEST) GO TO 320
      NSRCH(NCNT) = LN
   20 CONTINUE
   22 CONTINUE
      IF (NCNT.EQ.0) GO TO 200
      KCNT      = NCNT
      DO 25 K=1,KCNT
      KSRCH(K)  = NSRCH(K)
   25 CONTINUE
      GO TO 10
!
!     TETRAHEDRON FOUND. RE-INITIALIZE ARRAY NFLAG
!
   30 IF (IMON.EQ.0) GO TO 40
      DO 35 K=1,IMON
      NFLAG(IFLAG(K)) = 0
   35 CONTINUE
!
!     FIND ADDRESSES N1 AND N2 OF REMAINING TWO VERTICES OF TETRAHEDRON
!
   40 DO 45 K=1,4
      M         = NDC(K,LC)
      IF (M.NE.NA.AND.M.NE.NB) GO TO 50
   45 CONTINUE
   50 N1        = M
      N2        = NDC(1,LC)  +NDC(2,LC)  +NDC(3,LC)  +NDC(4,LC) &
                -NA  -NB  -N1
      IF (N1.GT.0) THEN
         IF (NVCNT(N1).EQ.4) THEN
            NVERT = N1
            RETURN
         ENDIF
      ENDIF
!
!     STARTING FROM FACE (N1,NA,NB) FIND SEQUENCE OF TETRAHEDRA
!     INCIDENT TO EDGE (NA,NB) UNTIL FACE (N2,NA,NB) IS REACHED
!
      NRBDY     = 0
      NRING     = 1

      CALL NEIGHB (LC,N1,NA,NB,JNEXT,K1,K2,NDC,NBH)

      NPNEXT    = NDC(1,JNEXT)  +NDC(2,JNEXT)  +NDC(3,JNEXT) &
                 +NDC(4,JNEXT)  -N1  -NA  -NB
      IRING(1)  = N1
      NTETKP(1) = JNEXT
      NFLAG(JNEXT) = 1
      IFLAG(1)  = JNEXT
      IF (IPROT(JNEXT).EQ.1) NRBDY = 1
   60 JPRE      = JNEXT
      NPRE      = NPNEXT

      CALL NEIGHB (JPRE,NPRE,NA,NB,JNEXT,K1,K2,NDC,NBH)

      NPNEXT    = NDC(1,JNEXT)  +NDC(2,JNEXT)  +NDC(3,JNEXT) &
                 +NDC(4,JNEXT)  -NPRE  -NA  -NB
      IF (NPRE.GT.0) THEN
         IF (NVCNT(NPRE).EQ.4) THEN
            NVERT = NPRE
            DO 62 K=1,NRING
            NFLAG(IFLAG(K)) = 0
   62       CONTINUE
            RETURN
         ENDIF
      ENDIF
      NRING     = NRING +1
      IF (NRING.GT.MXRING) GO TO 340
      IRING(NRING) = NPRE
      NTETKP(NRING) = JNEXT
      NFLAG(JNEXT) = 1
      IFLAG(NRING) = JNEXT
      IF (IPROT(JNEXT).EQ.1) NRBDY = 1
      IF (NPNEXT.EQ.N1) GO TO 350
      IF (NPNEXT.NE.N2) GO TO 60
      IF (N2.GT.0) THEN
         IF (NVCNT(N2).EQ.4) THEN
            NVERT = N2
            DO 64 K=1,NRING
            NFLAG(IFLAG(K)) = 0
   64       CONTINUE
            RETURN
         ENDIF
      ENDIF
      NRING     = NRING +1
      IF (NRING.GT.MXRING) GO TO 340
      IRING(NRING) = N2
      NTETKP(NRING) = LC
      NFLAG(LC) = 1
      IFLAG(NRING) = LC
      IMON      = NRING
!
!     FIND ALL TETRAHEDRA THAT ARE INCIDENT TO THE VERTICES NA AND NB
!
      NABDY     = 0
      NBBDY     = 0
      NP        = NA
   65 KSRCH(1)  = NPTET(NP)
      KCNT      = 1
      LC        = KSRCH(1)
      N1        = NDC(1,LC)
      N2        = NDC(2,LC)
      N3        = NDC(3,LC)
      N4        = NDC(4,LC)
      IF (LC.LE.0) GO TO 300
      IF (NFLAG(LC).GT.0) GO TO 70
      IF (N1.NE.NP.AND.N2.NE.NP.AND.N3.NE.NP.AND.N4.NE.NP) GO TO 70
      NFLAG(LC) = NP
      IMON      = IMON  +1
      IF (IMON.GT.MXNODE) GO TO 310
      IFLAG(IMON) = LC
      IF (IPROT(LC).EQ.1.AND.NP.EQ.NA) NABDY = 1
      IF (IPROT(LC).EQ.1.AND.NP.EQ.NB) NBBDY = 1
   70 NCNT      = 0
      DO 75 K=1,KCNT
      L         = KSRCH(K)
      DO 75 M=1,4
      LC        = NBH(M,L)
      IF (LC.LE.0) GO TO 300
      IF (NFLAG(LC).GT.0) GO TO 75
      N1        = NDC(1,LC)
      N2        = NDC(2,LC)
      N3        = NDC(3,LC)
      N4        = NDC(4,LC)
      IF (N1.NE.NP.AND.N2.NE.NP.AND.N3.NE.NP.AND.N4.NE.NP) GO TO 75
      NCNT      = NCNT  +1
      IF (NCNT.GT.MXTEST) GO TO 320
      NSRCH(NCNT) = LC
      NFLAG(LC) = NP
      IMON      = IMON  +1
      IF (IMON.GT.MXNODE) GO TO 310
      IFLAG(IMON) = LC
      IF (IPROT(LC).EQ.1.AND.NP.EQ.NA) NABDY = 1
      IF (IPROT(LC).EQ.1.AND.NP.EQ.NB) NBBDY = 1
   75 CONTINUE
      IF (NCNT.EQ.0) GO TO 85
      KCNT      = NCNT
      DO 80 K=1,KCNT
      KSRCH(K)  = NSRCH(K)
   80 CONTINUE
      GO TO 70
   85 IF (NP.EQ.NB) GO TO 90
      NP        = NB
      GO TO 65
!
!     FIND NEIGHBORS OF THE TETRAHEDRAL RING
!
   90 RADMAX    = 0.
      RATMAX    = 0.
      DO 115 I=1,NRING
      RADMAX    = MAX(RADMAX,RC(NTETKP(I)))
      RATMAX    = MAX(RATMAX,RAT(NTETKP(I)))
      IP1       = MOD(I,NRING)  +1

      CALL NEIGHB (NTETKP(I),IRING(I),IRING(IP1),NA,LNBR(I),KK,ISHK(I), &
                   NDC,NBH)

      CALL NEIGHB (NTETKP(I),IRING(I),IRING(IP1),NB,MNBR(I),KK,KSHK(I), &
                   NDC,NBH)

  115 CONTINUE
!
!     DEFINE COORDINATES OF NEW POINT TO REPLACE POINTS AT NA AND NB
!
!     Exit from COLAPS if edge (NA,NB) lies on the boundary.
!
      IF (NABDY.EQ.1.AND.NBBDY.EQ.1) GO TO 230

      NPASS   = 0
      XPT     = .5*(X(1,NA)  +X(1,NB))
      YPT     = .5*(X(2,NA)  +X(2,NB))
      ZPT     = .5*(X(3,NA)  +X(3,NB))
      ITYPT   = MAX(ITYP(NA),ITYP(NB))
      IF (NRBDY.EQ.NABDY.AND.NRBDY.EQ.NBBDY) GO TO 118
      XPT     = X(1,NA)
      YPT     = X(2,NA)
      ZPT     = X(3,NA)
      ITYPT   = ITYP(NA)
      IF (NRBDY.EQ.0.AND.NABDY.EQ.1.AND.NBBDY.EQ.0) GO TO 118
      XPT     = X(1,NB)
      YPT     = X(2,NB)
      ZPT     = X(3,NB)
      ITYPT   = ITYP(NB)
      NEXCH   = NA
      NA      = NB
      NB      = NEXCH
      IF (NRBDY.EQ.0.AND.NABDY.EQ.0.AND.NBBDY.EQ.1) THEN
         NABDY    = 1
         NBBDY    = 0
         GO TO 118
      ELSE
         GO TO 230
      ENDIF
!
!     CHECK WHETHER NEW TETRAHEDRAL RING FORMS A CONVEX ENSEMBLE
!
  118 NRNGP1  = NRING  +1
      DO 125 I=NRNGP1,IMON
      L       = IFLAG(I)
      IF (IPROT(L).EQ.1) GO TO 125
      N1      = NDC(1,L)
      N2      = NDC(2,L)
      N3      = NDC(3,L)
      N4      = NDC(4,L)
  120 IF (N4.EQ.NA.OR.N4.EQ.NB) GO TO 123
      NKEEP     = N1
      N1        = N2
      N2        = N3
      N3        = N4
      N4        = NKEEP
      GO TO 120
  123 RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      VTET1   = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                         +RNZ*(X(3,N4)  -X(3,N1))
      VTET2   = RNX*(XPT  -X(1,N1))  +RNY*(YPT  -X(2,N1)) &
                                     +RNZ*(ZPT  -X(3,N1))
      IF (ABS(VTET2).LT.TOLS.OR.VTET1*VTET2.LT.0.0) THEN
         IF (NPASS.EQ.1) THEN
             GO TO 220
         ELSE
             NPASS    = 1
             IF (N4.EQ.NA.AND.NABDY.EQ.0.AND.NBBDY.EQ.1) GO TO 220
             IF (N4.EQ.NB.AND.NABDY.EQ.1.AND.NBBDY.EQ.0) GO TO 220
             XPT      = X(1,N4)
             YPT      = X(2,N4)
             ZPT      = X(3,N4)
             ITYPT    = ITYP(N4)
             GO TO 118
         ENDIF
      ENDIF

      CALL CIRCUM2 (X,N1,N2,N3,XPT,YPT,ZPT,  &
                    XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (RADC.GT.50.*RADMAX) GO TO 220
      AREA      = TETAR3 (N1,N2,N3,XPT,YPT,ZPT,X)
      RCRIN     = RADC*AREA/VL
      IF (RCRIN.GT.50.*RATMAX) GO TO 220
  125 CONTINUE
!
!     MODIFY BOUNDARY DATA STRUCTURE IF THE EDGE (NA,NB) IS
!     A BOUNDARY EDGE
!
!
!     ASSIGN NEW POINT TO ADDRESS NA
!
!
      CALL OCTRMV (NA,X,NOCTR,NLINK,IDONE,NREF)

      CALL OCTRMV (NB,X,NOCTR,NLINK,IDONE,NREF)

      ITYP(NB) = -1
      X(1,NA) = XPT
      X(2,NA) = YPT
      X(3,NA) = ZPT
      ITYP(NA) = ITYPT

      CALL OCTFIL (NA,X,NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR)

      IDONE(NA) = 1
      NCOL    = NCOL  +1
!
!     REMOVE RING OF TETRAHEDRA AND UPDATE NBH ARRAY
!
      DO 130 I=1,NRING
      NBH(1,NTETKP(I)) = 0
      NBH(2,NTETKP(I)) = 0
      NBH(3,NTETKP(I)) = 0
      NBH(4,NTETKP(I)) = 0
      NBH(ISHK(I),LNBR(I)) = MNBR(I)
      NBH(KSHK(I),MNBR(I)) = LNBR(I)
      IF (IRING(I).GT.0) NVCNT(IRING(I)) = NVCNT(IRING(I))  -1
  130 CONTINUE
      NVCNT(NA)        = NVCNT(NA)  +NVCNT(NB)  -NRING  -2
!
!     REMOVE EDGE (NA,NB) AND EDGES JOINING NB TO RING POINTS
!     FROM DATA STRUCTURE NDG
!
      L1      = MIN(NA,NB)
      L2      = MAX(NA,NB)
      IEDG    = IDGP(L1)
  135 IF (IEDG.EQ.0) GO TO 330
      IF (L2.EQ.NDG(2,IEDG)) GO TO 140
      IEDG    = NDGP(IEDG)
      GO TO 135

  140 CALL EDGERM (IEDG,NDG,IDGP,NDGP)

      DO 152 I=1,NRING
      IF (IRING(I).LT.0) GO TO 152
      L1      = MIN(NB,IRING(I))
      L2      = MAX(NB,IRING(I))
      IEDG    = IDGP(L1)
  145 IF (IEDG.EQ.0) GO TO 330
      IF (L2.EQ.NDG(2,IEDG)) GO TO 150
      IEDG    = NDGP(IEDG)
      GO TO 145
!
  150 CALL EDGERM (IEDG,NDG,IDGP,NDGP)
!
  152 CONTINUE
!
!     UPDATE NDC ARRAYS WHICH REFER TO POINT NB
!
      DO 162 I=NRNGP1,IMON
      L       = IFLAG(I)
      IF (NFLAG(L).EQ.NA) GO TO 162
      N1      = NDC(1,L)
      N2      = NDC(2,L)
      N3      = NDC(3,L)
      N4      = NDC(4,L)
      IF (N1.EQ.NB) NDC(1,L) = NA
      IF (N2.EQ.NB) NDC(2,L) = NA
      IF (N3.EQ.NB) NDC(3,L) = NA
      IF (N4.EQ.NB) NDC(4,L) = NA
  162 CONTINUE
!
!     UPDATE NDG ARRAYS WHICH REFER TO POINT NB
!
      DO 177 I=NRNGP1,IMON
      L       = IFLAG(I)
      IF (NFLAG(L).EQ.NA) GO TO 177
      DO 175 K=1,4
      N       = NDC(K,L)
      IF (N.LT.0.OR.N.EQ.NA) GO TO 175
      L1      = MIN(N,NB)
      L2      = MAX(N,NB)
      IEDG    = IDGP(L1)
  165 IF (IEDG.EQ.0) GO TO 175
      IF (L2.EQ.NDG(2,IEDG)) GO TO 170
      IEDG    = NDGP(IEDG)
      GO TO 165

  170 CALL EDGERM (IEDG,NDG,IDGP,NDGP)

      M1      = MIN(N,NA)
      M2      = MAX(N,NA)
      INAD    = IDGP(M1)
      IF (INAD.EQ.0) GO TO 174
  172 IF (NDGP(INAD).EQ.0) GO TO 174
      INAD    = NDGP(INAD)
      GO TO 172
  174 NDG(1,IEDG) = M1
      NDG(2,IEDG) = M2
      NDGP(IEDG) = 0
      IF (INAD.EQ.0) IDGP(M1) = IEDG
      IF (INAD.GT.0) NDGP(INAD) = IEDG
  175 CONTINUE
  177 CONTINUE
!
!     COMPUTE VOLUME,CIRCUMCENTER AND CIRCUMRADIUS FOR NEW CELLS
!
      DO 195 K=NRNGP1,IMON
      L       = IFLAG(K)
      IF (IPROT(L).EQ.1) GO TO 195
      N1      = NDC(1,L)
      N2      = NDC(2,L)
      N3      = NDC(3,L)
      N4      = NDC(4,L)

      CALL CIRCUM (X,N1,N2,N3,N4,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (ISMALL.EQ.1) GO TO 360
      XCEN(L) = XCN
      YCEN(L) = YCN
      ZCEN(L) = ZCN
      VOL(L)  = VL
      RC(L)   = RADC
      AREA    = TETAR (L,X,NDC)
      RAT(L)  = RC(L)*AREA/VOL(L)
      DO 190 I=1,4
      N        = NDC(I,L)
      NPTET(N) = L
  190 CONTINUE
  195 CONTINUE
      GO TO 250
!
!     FAILURE MODES
!
  200 NFAIL    = 1
!     WRITE (6,700) NA,NB
  700 FORMAT('EDGE WITH VERTICES ',I6,' , ',I6,' DOES NOT EXIST.'/ &
             'COLLAPSE OF THESE TWO POINTS IS NOT POSSIBLE.')
      GO TO 250
  220 NFAIL    = 3
!     WRITE (6,720) NA,NB
  720 FORMAT('ATTEMPTED COLLAPSE OF EDGE WITH VERTICES ',I6,' , ',I6/ &
             'CREATES A NON-CONVEX TETRAHEDRON')
      GO TO 250
  230 NFAIL    = 4
!     WRITE (6,725) NA,NB
  725 FORMAT('ATTEMPTED COLLAPSE OF AN EDGE THAT DOES NOT LIE IN THE'/ &
             'BOUNDARY SURFACE BUT WHOSE END-POINTS ',I6,' AND ',I6/ &
             'BOTH LIE IN THE BOUNDARY SURFACE')
!
!     RE-INITIALIZE NFLAG ARRAY
!
  250 DO 260 I=1,IMON
      NFLAG(IFLAG(I)) = 0
  260 CONTINUE
      RETURN
  300 WRITE (6,600) NA,NB,LC
      STOP
  310 WRITE (6,610)
      STOP
  320 WRITE (6,620)
      STOP
  330 WRITE (6,630)
      STOP
  340 WRITE (6,640)
      STOP
  350 WRITE (6,650)
      STOP
  360 WRITE (6,660)
      STOP
  600 FORMAT(///5X,'INVALID TETRAHEDRON ADDRESS FOUND WHILE SEARCHING'/ &
                5X,'FOR A TETRAHEDRON INCIDENT TO EDGE WITH VERTICES', &
                   I6,' , ',I6/5X,'TETRAHEDRON ADDRESS IS ',I7/ &
                5X,'PROGRAM STOPPED IN ROUTINE COLAPS')
  610 FORMAT(///5X,'DIMENSION OF ARRAY IFLAG EXCEEDED. INCREASE SIZE'/ &
                5X,'OF MNODE. PROGRAM STOPPED IN COLAPS')
  620 FORMAT(///5X,'DIMENSION OF ARRAY NSRCH EXCEEDED. INCREASE SIZE'/ &
                5X,'OF MTEST. PROGRAM STOPPED IN COLAPS')
  650 FORMAT(///5X,'SEARCH FOR TETRAHEDRAL RING HAS RETURNED TO THE'/ &
                5X,'STARTING FACE. THIS INDICATES AN INCONSISTENCY IN'/ &
                5X,'THE TETRAHEDRAL ENSEMBLE.'/ &
                5X,'PROGRAM STOPPED IN COLAPS')
  630 FORMAT(///5X,'UNABLE TO FIND EDGE ADDRESS FOR A NEW TETRAHEDRON'/ &
                5X,'PROGRAM STOPPED IN ROUTINE COLAPS') 
  640 FORMAT(///5X,'DIMENSION OF ARRAY IRING EXCEEDED. INCREASE SIZE'/ &
                5X,'OF MRING. PROGRAM STOPPED IN COLAPS')
  660 FORMAT(///5X,'AT LEAST ONE NEW TETRAHEDRON HAS TOO SMALL A VOLUME' &
               /5X,'PROGRAM STOPPED IN COLAPS')
      END SUBROUTINE COLAPS






!
!     ******************************************************************
!
      SUBROUTINE SMOOTH (X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                         NDG,IDGP,NDGP,NEDGE, &
                         VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                         IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                         XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                         XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                         IPOINT,NPOINT,IFLAG,NFLAG, &
                         DX,DY,DZ,DS,VLT,IRING,NTETKP,NFAD,NEWC, &
                         NBHKP,IEDKP,LNBR,ISHK,MNBR,KSHK,NPP, &
                         NFILL,NEWCEL,NTRI, &
                         NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                         NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                         LISTF,VOLMIN,RCMX,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  OPTIMIZE MESH QUALITY BY ALTERNATELY SWAPPING EDGE/FACE       *
!     *  COMBINATIONS AND ADJUSTING POINT POSITIONS.                   *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,NCELL,NEDGE,NNODE
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: IDONE(*),NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: ITYP(*),NPTET(*),NACPT(*),LISTF(*)
      INTEGER :: NPOINT(*),IPOINT(*)
      INTEGER :: IFLAG(*),NFLAG(*)
      INTEGER :: LNKUP(*),LNKDN(*)
      INTEGER :: NSRCH(*),KSRCH(*)
      INTEGER :: IRING(*),NTETKP(*),NFAD(3,*),NEWC(*), &
                 NBHKP(3,*),IEDKP(4,*),LNBR(*),ISHK(*), &
                 MNBR(*),KSHK(*),NPP(*)
      INTEGER :: NTRI(3,*),NFILL(*),NEWNBH(4,*),NOLD(*), &
                 NEWCEL(*),NSHAKE(*),NCAV(4,*)
      INTEGER :: NEDGRM(*)
      INTEGER :: LDEL(*),NCAVFC(3,*),IKEEP(*)
      DOUBLE PRECISION :: VOLMIN,RCMX,TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XOCTR(2,*),YOCTR(2,*),ZOCTR(2,*), &
                          XHOLD(2,*),YHOLD(2,*),ZHOLD(2,*), &
                          XFAR(2),YFAR(2),ZFAR(2),XKEEP(2), &
                          YKEEP(2),ZKEEP(2)
      DOUBLE PRECISION :: DX(*),DY(*),DZ(*),DS(*),VLT(*)
      DOUBLE PRECISION :: XC(*),YC(*),ZC(*),V(*),RAD(*),RCRIN(*)
      DOUBLE PRECISION :: ANG(6)
      
      INTEGER :: I,I1,I2,I3,J,K,L,LOOP,M1,M2,N,NA,NB,NCHK,NFAIL, &
                 NMAX,NMIN,NN,NP,NPASS,NPT,NRING,NSUM,NTET, &
                 NWAIT,N1,N2,N3,N4
      DOUBLE PRECISION :: ANGMAX,ANGMX1,ANGMX2,ANGMX3
      
!
!     ******************************************************************
!
!     CREATE ORDERED LIST OF CELLS ACCORDING TO CIRCUM-RADIUS TO
!     IN-RADIUS RATIO
!
      WRITE (6,610)

      DO 140 LOOP=1,10

      CALL TREE (RAT,LNKDN,LNKUP,NFLAG,LISTF,NACPT,NWAIT, &
                 NBH,IPROT,NCELL)

   30 WRITE (6,910) NWAIT,LOOP
  910 FORMAT(//'IN SMOOTH,  NWAIT = ',I7,'   ITERATION ',I2)
!
!     SEARCH THROUGH LIST OF CELLS AND TRY TO REMOVE SLIVERS BY
!     EDGE/FACE SWAPS
!
      NCHK   = 0
      L      = 1
   40 J      = LISTF(L)
      IF (NBH(1,J).LE.0) GO TO 130
      IF (NACPT(J).EQ.1) GO TO 130
      NACPT(J) = 1
      N1     = NDC(1,J)
      N2     = NDC(2,J)
      N3     = NDC(3,J)
      N4     = NDC(4,J)
      ANG(1)   = DIHED (N1,N2,N3,N4,X)
      IF (ANG(1).GT.180.) ANG(1) = 360.  -ANG(1)
      ANG(2)   = DIHED (N2,N3,N4,N1,X)
      IF (ANG(2).GT.180.) ANG(2) = 360.  -ANG(2)
      ANG(3)   = DIHED (N3,N4,N1,N2,X)
      IF (ANG(3).GT.180.) ANG(3) = 360.  -ANG(3)
      ANG(4)   = DIHED (N4,N1,N2,N3,X)
      IF (ANG(4).GT.180.) ANG(4) = 360.  -ANG(4)
      ANG(5)   = DIHED (N1,N3,N2,N4,X)
      IF (ANG(5).GT.180.) ANG(5) = 360.  -ANG(5)
      ANG(6)   = DIHED (N2,N4,N1,N3,X)
      IF (ANG(6).GT.180.) ANG(6) = 360.  -ANG(6)
!
!     FIND THE THREE EDGES WITH THE LARGEST DIHEDRAL ANGLES
!
      I1     = 1
      ANGMX1 = ANG(1)
      DO 50 K=2,6
      IF (ANG(K).LT.ANGMX1) GO TO 50
      I1     = K
      ANGMX1 = ANG(K)
   50 CONTINUE
      I2     = 0
      DO 60 K=1,6
      IF (K.EQ.I1) GO TO 60
      IF (I2.EQ.0) GO TO 55
      IF (ANG(K).LT.ANGMX2) GO TO 60
   55 I2     = K
      ANGMX2 = ANG(K)
   60 CONTINUE
      I3     = 0
      DO 70 K=1,6
      IF (K.EQ.I1.OR.K.EQ.I2) GO TO 70
      IF (I3.EQ.0) GO TO 65
      IF (ANG(K).LT.ANGMX3) GO TO 70
   65 I3     = K
      ANGMX3 = ANG(K)
   70 CONTINUE
!
!     DETERMINE WHETHER THE SINGULAR TETRAHEDRON IS A FLAT OR A SLIVER
!
      NP     = 0
      NMIN   = MIN(I1,I2,I3)
      NMAX   = MAX(I1,I2,I3)
      NSUM   = I1  +I2  +I3
      IF (NMIN.EQ.1.AND.NMAX.EQ.5.AND.NSUM.EQ.10) NP = N1
      IF (NMIN.EQ.1.AND.NMAX.EQ.6.AND.NSUM.EQ.9)  NP = N2
      IF (NMIN.EQ.2.AND.NMAX.EQ.5.AND.NSUM.EQ.10) NP = N3
      IF (NMIN.EQ.3.AND.NMAX.EQ.6.AND.NSUM.EQ.13) NP = N4
      IF (NP.EQ.0.OR.ANGMX3.LT.120.) GO TO 110
!
!     A FLAT TETRAHEDRON HAS BEEN FOUND
!
   75 IF (NP.EQ.N4) GO TO 80
      NN     = N1
      N1     = N2
      N2     = N3
      N3     = N4
      N4     = NN
      GO TO 75

   80 CALL TRISWP (N4,N1,N2,N3,J,NFAIL, &
                   X,NDC,NBH,IPROT,NCELL,NDG,IDGP,NDGP,NEDGE, &
                   VOL,XCEN,YCEN,ZCEN,RC,RAT,NPTET,NACPT,TOLV)

      IF (NFAIL.EQ.0) THEN
         GO TO 125
      ELSE
         GO TO 128
      ENDIF
!
!     A SLIVER HAS BEEN FOUND
!
  110 IF (MAX(ANGMX1,ANGMX2).LT.120.) GO TO 135
      NPASS  = 0
      I      = I1
      ANGMAX = ANGMX1
  115 IF (ANGMAX.LT.120.) GO TO 130
      NA     = N1
      NB     = N2
      M1     = N3
      IF (I.EQ.1) GO TO 120
      NA     = N2
      NB     = N3
      M1     = N4
      IF (I.EQ.2) GO TO 120
      NA     = N3
      NB     = N4
      M1     = N1
      IF (I.EQ.3) GO TO 120
      NA     = N4
      NB     = N1
      M1     = N2
      IF (I.EQ.4) GO TO 120
      NA     = N1
      NB     = N3
      M1     = N4
      IF (I.EQ.5) GO TO 120
      NA     = N2
      NB     = N4
      M1     = N1
  120 M2     = N1  +N2  +N3  +N4  -NA  -NB  -M1

      CALL REPLACE (M1,M2,NA,NB,J,NRING,NPASS,NFAIL, &
                    X,DENS,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                    NDG,IDGP,NDGP,NEDGE,IPOINT,NPOINT, &
                    VOL,XCEN,YCEN,ZCEN,RC,RAT,NPTET,NACPT, &
                    NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR, &
                    IRING,NTETKP,LNBR,ISHK,MNBR,KSHK, &
                    NFAD,NBHKP,IEDKP,NEWC,DX,DY,DZ,DS,NPP,VLT,TOLV)

      IF (NFAIL.EQ.0) GO TO 125
      IF (NPASS.EQ.1.OR.NFAIL.EQ.2) GO TO 128
      NPASS    = 1
      I        = I2
      ANGMAX   = ANGMX2
      GO TO 115
  125 NCHK     = NCHK  +1
      GO TO 40
  128 IF (LOOP.LE.4) THEN

         CALL PUTPNT (J,NFAIL, &
                      X,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                      NDG,IDGP,NDGP,NEDGE, &
                      VOL,XCEN,YCEN,ZCEN,RC,RAT,DENS,NPTET,NACPT, &
                      IDONE,NREF,NLINK,NOCTR,IOCTR,XFAR,YFAR,ZFAR, &
                      XOCTR,YOCTR,ZOCTR,XHOLD,YHOLD,ZHOLD, &
                      XKEEP,YKEEP,ZKEEP,KSRCH,NSRCH, &
                      IPOINT,NPOINT,IFLAG,NFLAG,NFILL,NEWCEL,NTRI, &
                      NCAV,NSHAKE,NEWNBH,NOLD,NCAVFC,IKEEP,LDEL, &
                      NEDGRM,XC,YC,ZC,V,RAD,RCRIN,LNKUP,LNKDN, &
                      VOLMIN,RCMX,TOLV)

         IF (NFAIL.EQ.0) GO TO 125
      ENDIF
  130 L      = L  +1
      IF (L.GT.NWAIT) GO TO 135
      GO TO 40
  135 CONTINUE
      WRITE (6,950) NCHK
  950 FORMAT(I6,' SLIVERS REMOVED')
!
!     CYCLE THROUGH FIELD POINTS AND TRY ADJUST POINT POSITIONS
!     BY A LAPLACIAN SMOOTHER
!
      IF (NCHK.EQ.0) GO TO 145
  140 CONTINUE
  145 NTET     = 0
      DO 160 L=1,NCELL
      IF (IPROT(L).EQ.1) GO TO 160
      IF (NBH(1,L).EQ.0) GO TO 160
      NTET     = NTET  +1
  160 CONTINUE
      NPT      = 0
      DO 165 N=1,NNODE
      IF (ITYP(N).LT.0) GO TO 165
      NPT      = NPT  +1
  165 CONTINUE
      WRITE (6,600) NPT,NTET
      RETURN
  600 FORMAT(/'MESH OPTIMIZATION COMPLETE'/ &
            5X,I7,' MESH POINTS'/ &
            5X,I7,' MESH CELLS'/)
  610 FORMAT(/'BEGIN MESH OPTIMIZATION')
      END SUBROUTINE SMOOTH





!
!     ******************************************************************
! 
      SUBROUTINE REPLACE (N1,N2,NA,NB,J,NRING,NPASS,NFAIL, &
                          X,DENS,ITYP,NNODE,NDC,NBH,IPROT,NCELL, &
                          NDG,IDGP,NDGP,NEDGE,IPOINT,NPOINT, &
                          VOL,XCEN,YCEN,ZCEN,RC,RAT,NPTET,NACPT, &
                          NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR, &
                          IRING,NTETKP,LNBR,ISHK,MNBR,KSHK, &
                          NFAD,NBHKP,IEDKP,NEWC,DX,DY,DZ,DS,NP,VLT,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  FIND THE RING OF TETRAHEDRA SURROUNDING EDGE (NA,NB) AND      *
!     *  CREATE NEW ENSEMBLE OF TETRAHEDRA IN ORDER TO IMPROVE MESH    *
!     *  QUALITY                                                       *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,J,NA,NB,NCELL,NEDGE,NFAIL,NNODE,NPASS,NRING,N1,N2
      INTEGER :: ITYP(*),NPOINT(*),IPOINT(*),NPTET(*),NACPT(*)
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*),IDGP(*),NDGP(*)
      INTEGER :: NREF(*),NLINK(*),NOCTR(2,*)
      INTEGER :: IRING(*),NTETKP(*),NFAD(3,*),NEWC(*), &
                 NBHKP(3,*),IEDKP(4,*),LNBR(*),ISHK(*), &
                 MNBR(*),KSHK(*),NP(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: VOL(*),XCEN(*),YCEN(*),ZCEN(*),RC(*),RAT(*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2)
      DOUBLE PRECISION :: DX(*),DY(*),DZ(*),DS(*),VLT(*)

      INTEGER :: I,IEDG,IP1,ISMALL,JNEXT,JPRE,J1,J2,K,KCNT,KEND,KHALF, &
                 KK,KMAX,KMIN,KM1,KPLUS,KPLUSM,KPLUSP,KP1,KSUM,K1,K2, &
                 K3,K4,L,LOOP,L1,L2,M,MEND,MMAX,MMIN,MSUM,M1,M2,M3,M4, &
                 M5,N,NEDK,NPNEXT,NPRE,NRNGP1
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,RADC,RADMAX,RATMAX, &
                          RCRIN,RNX,RNY,RNZ,VL,VMIN,VORIG,VTET1,VTET2, &
                          XCN,YCN,ZCN
!
!     ******************************************************************
!
      NFAIL     = 0
      VORIG     = VOL(J)
!
!     STARTING FROM SLIVER FACE (N1,NA,NB) FIND SEQUENCE OF TETRAHEDRA
!     INCIDENT TO EDGE (NA,NB) UNTIL FACE (N2,NA,NB) IS REACHED
!
      NRING     = 1

      CALL NEIGHB (J,N1,NA,NB,JNEXT,K1,K2,NDC,NBH)

      NPNEXT    = NDC(1,JNEXT)  +NDC(2,JNEXT)  +NDC(3,JNEXT) &
                 +NDC(4,JNEXT)  -N1  -NA  -NB
      IRING(1)  = N1
      NTETKP(1) = JNEXT
   10 JPRE      = JNEXT
      NPRE      = NPNEXT

      CALL NEIGHB (JPRE,NPRE,NA,NB,JNEXT,K1,K2,NDC,NBH)

      NPNEXT    = NDC(1,JNEXT)  +NDC(2,JNEXT)  +NDC(3,JNEXT) &
                 +NDC(4,JNEXT)  -NPRE  -NA  -NB
      NRING     = NRING +1
      IF (NRING.GT.MXRING) GO TO 310
      IRING(NRING) = NPRE
      NTETKP(NRING) = JNEXT
      IF (NPNEXT.EQ.N1) GO TO 300
      IF (NPNEXT.NE.N2) GO TO 10
      NRING     = NRING +1
      IRING(NRING) = N2
      NTETKP(NRING) = J
!
!     FIND NEIGHBORS OF THE TETRAHEDRAL RING
!
      RADMAX    = 0.
      RATMAX    = 0.
      DO 20 I=1,NRING
      IF (IPROT(NTETKP(I)).EQ.1) GO TO 200
      RADMAX    = MAX(RADMAX,RC(NTETKP(I)))
      RATMAX    = MAX(RATMAX,RAT(NTETKP(I)))
      IP1        = MOD(I,NRING)  +1

      CALL NEIGHB (NTETKP(I),IRING(I),IRING(IP1),NA,LNBR(I),KK,ISHK(I), &
                   NDC,NBH)
                                                         
      CALL NEIGHB (NTETKP(I),IRING(I),IRING(IP1),NB,MNBR(I),KK,KSHK(I), &
                   NDC,NBH)

   20 CONTINUE
      KHALF       = NRING  -2

      CALL TESSEL (IRING,NRING,NFAD,NBHKP,IEDKP,NEDK,0,0, &
                   X,IPOINT,NPOINT,DX,DY,DZ,DS,NP)
!
!     CHECK WHETHER TETRAHEDRAL RING FORMS A CONVEX ENSEMBLE
!
      DO 30 K=1,KHALF
      M1      = NFAD(1,K)
      M2      = NFAD(2,K)
      M3      = NFAD(3,K)

      CALL CIRCUM (X,M1,M2,M3,NA,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (RADC.GT.RADMAX) GO TO 200
      AREA      = TETAR2 (M1,M2,M3,NA,X)
      RCRIN     = RADC*AREA/VL
      IF (RCRIN.GT.RATMAX) GO TO 200

      CALL CIRCUM (X,M1,M2,M3,NB,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (RADC.GT.RADMAX) GO TO 200
      AREA      = TETAR2 (M1,M2,M3,NB,X)
      RCRIN     = RADC*AREA/VL
      IF (RCRIN.GT.RATMAX) GO TO 200
      RNX     = COFACT (X(2,M1),X(2,M2),X(2,M3),X(3,M1),X(3,M2),X(3,M3))
      RNY     = COFACT (X(3,M1),X(3,M2),X(3,M3),X(1,M1),X(1,M2),X(1,M3))
      RNZ     = COFACT (X(1,M1),X(1,M2),X(1,M3),X(2,M1),X(2,M2),X(2,M3))
      VTET1   = RNX*(X(1,NA)  -X(1,M1))  +RNY*(X(2,NA)  -X(2,M1)) &
                                         +RNZ*(X(3,NA)  -X(3,M1))
      VTET2   = RNX*(X(1,NB)  -X(1,M1))  +RNY*(X(2,NB)  -X(2,M1)) &
                                         +RNZ*(X(3,NB)  -X(3,M1))
      IF (ABS(VTET1).LT.VORIG.OR.ABS(VTET2).LT.VORIG) GO TO 100
      IF (VTET1*VTET2.GE.0.0) GO TO 100
   30 CONTINUE
      DO 39 LOOP=1,2
      DO 38 K=1,KHALF
      DO 38 I=1,3
      L       = NBHKP(I,K)
      IF (L.LT.0) GO TO 38
      J2      = LNBR(L)
      IF (LOOP.EQ.2) J2 = MNBR(L)
      IF (IPROT(J2).EQ.1) GO TO 38
      K1      = NFAD(1,K)
      K2      = NFAD(2,K)
      K3      = NFAD(3,K)
      K4      = NA
      IF (LOOP.EQ.2) K4 = NB
      KEND    = K4
      M1      = NDC(1,J2)
      M2      = NDC(2,J2)
      M3      = NDC(3,J2)
      M4      = NDC(4,J2)
   31 MEND    = M4
      KMIN    = MIN(K1,K2,K3)
      KMAX    = MAX(K1,K2,K3)
      KSUM    = K1  +K2  +K3
   32 MMIN    = MIN(M1,M2,M3)
      MMAX    = MAX(M1,M2,M3)
      MSUM    = M1  +M2  +M3
      IF (MMIN.EQ.KMIN.AND.MMAX.EQ.KMAX.AND.MSUM.EQ.KSUM) GO TO 34
      IF (M1.EQ.MEND) GO TO 33
      M       = M1
      M1      = M2
      M2      = M3
      M3      = M4
      M4      = M
      GO TO 32
   33 IF (K1.EQ.KEND) GO TO 340
      KK      = K1
      K1      = K2
      K2      = K3
      K3      = K4
      K4      = KK
      GO TO 31
   34 M5      = K4
      RNX     = COFACT (X(2,M1),X(2,M2),X(2,M3),X(3,M1),X(3,M2),X(3,M3))
      RNY     = COFACT (X(3,M1),X(3,M2),X(3,M3),X(1,M1),X(1,M2),X(1,M3))
      RNZ     = COFACT (X(1,M1),X(1,M2),X(1,M3),X(2,M1),X(2,M2),X(2,M3))
      VTET1   = RNX*(X(1,M4)  -X(1,M1))  +RNY*(X(2,M4)  -X(2,M1)) &
                                         +RNZ*(X(3,M4)  -X(3,M1))
      VTET2   = RNX*(X(1,M5)  -X(1,M1))  +RNY*(X(2,M5)  -X(2,M1)) &
                                         +RNZ*(X(3,M5)  -X(3,M1))
      IF (VTET1*VTET2.GE.0.) GO TO 100
   38 CONTINUE
   39 CONTINUE
!
!     ASSIGN ADDRESSES OF NEW CELLS
!
      KCNT        = 2*(NRING  -2)
      DO 40 K=1,KCNT
      IF (K.GT.NRING) GO TO 35
      NEWC(K)     = NTETKP(K)
      NACPT(NEWC(K)) = 1
      GO TO 40
   35 NCELL       = NCELL  +1
      NEWC(K)     = NCELL
      NACPT(NEWC(K)) = 1
   40 CONTINUE
!
!     ASSIGN ADDRESSES OF NEW NDC ARRAYS
!
      DO 60 K=1,KHALF
      NDC(1,NEWC(K)) = NFAD(1,K)
      NDC(2,NEWC(K)) = NFAD(2,K)
      NDC(3,NEWC(K)) = NFAD(3,K)
      NDC(4,NEWC(K)) = NA
      KPLUS           = K  +KHALF
      NDC(1,NEWC(KPLUS)) = NFAD(1,K)
      NDC(2,NEWC(KPLUS)) = NFAD(2,K)
      NDC(3,NEWC(KPLUS)) = NFAD(3,K)
      NDC(4,NEWC(KPLUS)) = NB
   60 CONTINUE
!
!     UPDATE NBH ARRAY
!
      DO 75 K=1,NRING
      NBH(1,NTETKP(K)) = 0
   75 CONTINUE
      DO 90 K=1,KHALF
      KPLUS      = K  +KHALF
      DO 85 I=1,3
      L          = NBHKP(I,K)
      IF (L.GT.0) GO TO 80
      L          = -L
      NBH(I,NEWC(K)) = NEWC(L)
      NBH(I,NEWC(KPLUS)) = NEWC(L+KHALF)
      GO TO 85
   80 NBH(I,NEWC(K)) = LNBR(L)
      NBH(I,NEWC(KPLUS)) = MNBR(L)
      NBH(ISHK(L),LNBR(L)) = NEWC(K)
      NBH(KSHK(L),MNBR(L)) = NEWC(KPLUS)
   85 CONTINUE
      NBH(4,NEWC(K)) = NEWC(KPLUS)
      NBH(4,NEWC(KPLUS)) = NEWC(K)
   90 CONTINUE
!
!     UPDATE NDG ARRAY
!
      IF (NRING.EQ.NEDK) GO TO 145
      NRNGP1     = NRING  +1
      DO 96 I=NRNGP1,NEDK
      M1         = IEDKP(1,I)
      M2         = IEDKP(2,I)
      L1         = MIN(M1,M2)
      L2         = MAX(M1,M2)
      IEDG       = IDGP(L1)
      IF (IEDG.EQ.0) GO TO 94
   92 IF (NDGP(IEDG).EQ.0) GO TO 94
      IEDG       = NDGP(IEDG)
      GO TO 92
   94 NEDGE      = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 370
      NDG(1,NEDGE) = L1
      NDG(2,NEDGE) = L2
      NDGP(NEDGE)  = 0
      IF (IEDG.EQ.0) IDGP(L1) = NEDGE
      IF (IEDG.GT.0) NDGP(IEDG) = NEDGE
   96 CONTINUE
      GO TO 145
!
!     TETRAHEDRAL ENSEMBLE IS NON-CONVEX. INSERT POINT INSIDE
!     THE ENSEMBLE TO REMOVE SLIVER
!
  100 CALL NEIGHB (J,N1,N2,NA,J1,K1,K2,NDC,NBH)

      CALL NEIGHB (J,N1,N2,NB,J2,K1,K2,NDC,NBH)

      IF (IPROT(J1).EQ.0.AND.IPROT(J2).EQ.0.AND.NPASS.EQ.0) GO TO 210

      CALL SLIVER (IRING,NRING,NA,NB,NFAIL,VORIG, &
                   X,DENS,ITYP,IPOINT,NNODE, &
                   NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR,VLT)

      IF (NFAIL.NE.0) RETURN
!
!     ASSIGN ADDRESSES OF NEW CELLS
!
      KCNT        = 2*NRING
      DO 110 K=1,KCNT
      IF (K.GT.NRING) GO TO 105
      NEWC(K)     = NTETKP(K)
      NACPT(NEWC(K)) = 1
      GO TO 110
  105 NCELL       = NCELL  +1
      NEWC(K)     = NCELL
      NACPT(NEWC(K)) = 1
  110 CONTINUE
!
!     ASSIGN ADDRESSES OF NEW NDC ARRAYS
!
      DO 120 K=1,NRING
      KP1         = MOD(K,NRING)  +1
      M1          = IRING(K)
      M2          = IRING(KP1)
      NDC(1,NEWC(K)) = M1
      NDC(2,NEWC(K)) = M2
      NDC(3,NEWC(K)) = NNODE
      NDC(4,NEWC(K)) = NA
      IPROT(NEWC(K)) = 0
      KPLUS           = K  +NRING
      NDC(1,NEWC(KPLUS)) = M1
      NDC(2,NEWC(KPLUS)) = M2
      NDC(3,NEWC(KPLUS)) = NNODE
      NDC(4,NEWC(KPLUS)) = NB
      IPROT(NEWC(KPLUS)) = 0
  120 CONTINUE
!
!     UPDATE NBH ARRAY
!
      DO 125 K=1,NRING
      NBH(1,NTETKP(K)) = 0
  125 CONTINUE
      DO 130 K=1,NRING
      KM1        = MOD(NRING-2+K,NRING)  +1
      KP1        = MOD(K,NRING)  +1
      KPLUS      = K  +NRING
      KPLUSM     = KM1  +NRING
      KPLUSP     = KP1  +NRING
      NBH(1,NEWC(K)) = NEWC(KM1)
      NBH(2,NEWC(K)) = NEWC(KP1)
      NBH(3,NEWC(K)) = LNBR(K)
      NBH(4,NEWC(K)) = NEWC(KPLUS)
      NBH(ISHK(K),LNBR(K)) = NEWC(K)
      NBH(1,NEWC(KPLUS)) = NEWC(KPLUSM)
      NBH(2,NEWC(KPLUS)) = NEWC(KPLUSP)
      NBH(3,NEWC(KPLUS)) = MNBR(K)
      NBH(4,NEWC(KPLUS)) = NEWC(K)
      NBH(KSHK(K),MNBR(K)) = NEWC(KPLUS)
  130 CONTINUE
!
!     UPDATE NDG ARRAY
!
      DO 136 K=1,NRING
      IEDG       = IDGP(IRING(K))
      IF (IEDG.EQ.0) GO TO 134
  132 IF (NDGP(IEDG).EQ.0) GO TO 134
      IEDG       = NDGP(IEDG)
      GO TO 132
  134 NEDGE      = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 370
      NDG(1,NEDGE) = IRING(K)
      NDG(2,NEDGE) = NNODE
      NDGP(NEDGE)  = 0
      IF (IEDG.EQ.0) IDGP(IRING(K)) = NEDGE
      IF (IEDG.GT.0) NDGP(IEDG) = NEDGE
  136 CONTINUE
      IEDG       = IDGP(NA)
      IF (IEDG.EQ.0) GO TO 140
  138 IF (NDGP(IEDG).EQ.0) GO TO 140
      IEDG       = NDGP(IEDG)
      GO TO 138
  140 NEDGE      = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 370
      NDG(1,NEDGE) = NA
      NDG(2,NEDGE) = NNODE
      NDGP(NEDGE)  = 0
      IF (IEDG.EQ.0) IDGP(NA) = NEDGE
      IF (IEDG.GT.0) NDGP(IEDG) = NEDGE
      IEDG       = IDGP(NB)
      IF (IEDG.EQ.0) GO TO 144
  142 IF (NDGP(IEDG).EQ.0) GO TO 144
      IEDG       = NDGP(IEDG)
      GO TO 142
  144 NEDGE      = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 370
      NDG(1,NEDGE) = NB
      NDG(2,NEDGE) = NNODE
      NDGP(NEDGE)  = 0
      IF (IEDG.EQ.0) IDGP(NB) = NEDGE
      IF (IEDG.GT.0) NDGP(IEDG) = NEDGE
!
!     REMOVE OLD EDGE JOINING POINTS NA AND NB
!
  145 L1      = MIN(NA,NB)
      L2      = MAX(NA,NB)
      IEDG    = IDGP(L1)
  146 IF (L2.EQ.NDG(2,IEDG)) GO TO 148
      IEDG    = NDGP(IEDG)
      GO TO 146

  148 CALL EDGERM (IEDG,NDG,IDGP,NDGP)
!
!     COMPUTE VOLUME,CIRCUMCENTER AND CIRCUMRADIUS FOR NEW CELLS
!
      DO 150 K=1,KCNT
      M1        = NDC(1,NEWC(K))
      M2        = NDC(2,NEWC(K))
      M3        = NDC(3,NEWC(K))
      M4        = NDC(4,NEWC(K))

      CALL CIRCUM (X,M1,M2,M3,M4,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (K.EQ.1) VMIN = VL
      IF (K.GT.1) VMIN = MIN(VMIN,VL)

      IF (ISMALL.EQ.1) WRITE (6,960) K,VL,RADC, &
                       M1,X(1,M1),X(2,M1),X(3,M1), &
                       M2,X(1,M2),X(2,M2),X(3,M2), &
                       M3,X(1,M3),X(2,M3),X(3,M3), &
                       M4,X(1,M4),X(2,M4),X(3,M4) 
  960 FORMAT('RING ',I2,' VOLUME = ',E13.5,' RC = ',E13.5/ &
             'M1 ',I6,' X = ',F6.3,' Y = ',F6.3,' Z = ',F6.3/ &
             'M2 ',I6,' X = ',F6.3,' Y = ',F6.3,' Z = ',F6.3/ &
             'M3 ',I6,' X = ',F6.3,' Y = ',F6.3,' Z = ',F6.3/ &
             'M4 ',I6,' X = ',F6.3,' Y = ',F6.3,' Z = ',F6.3)

      IF (ISMALL.EQ.1) GO TO 320
      XCEN(NEWC(K)) = XCN
      YCEN(NEWC(K)) = YCN
      ZCEN(NEWC(K)) = ZCN
      VOL(NEWC(K))  = VL
      RC(NEWC(K))   = RADC
      AREA       = TETAR (NEWC(K),X,NDC)
      RAT(NEWC(K)) = RC(NEWC(K))*AREA/VOL(NEWC(K))
!     WRITE (6,961) K,VL,RADC,RAT(NEWC(K))
! 961 FORMAT('NEW CELLS,  K = ',I2,' VOL ',E13.5,
!    .                             ' RC ',E13.5,' RAT ',E13.5)
  150 CONTINUE

!     WRITE (6,990) VORIG,VMIN
! 990 FORMAT('VORIG = ',E13.5,' VMIN = ',E13.5)
!
!     UPDATE NPTET ARRAY
!
      DO 155 K=1,KCNT
      DO 155 I=1,4
      N           = NDC(I,NEWC(K))
      NPTET(N)   = NEWC(K)
  155 CONTINUE
      RETURN
  200 NFAIL       = 1
!     WRITE (6,920)
! 920 FORMAT('AT LEAST ONE OF THE RING TETRAHEDRA IS PROTECTED')
      RETURN
  210 NFAIL       = 1
!     WRITE (6,921)
! 921 FORMAT('NON BODY STICKING SLIVER. TRY OPPOSITE DIRECTION')
      RETURN
  300 WRITE (6,600)
      STOP
  310 WRITE (6,610)
      STOP
  320 WRITE (6,620)
      STOP
  340 WRITE (6,640)
      STOP
  360 WRITE (6,660)
      STOP
  370 WRITE (6,670)
      STOP
  600 FORMAT(///5X,'SEARCH FOR TETRAHEDRAL RING HAS RETURNED TO THE'/ &
                5X,'STARTING FACE. THIS INDICATES AN INCONSISTENCY IN'/ &
                5X,'THE TETRAHEDRAL ENSEMBLE.'/ &
                5X,'PROGRAM STOPPED IN REPLACE')
  610 FORMAT(///5X,'DIMENSION OF ARRAY IRING EXCEEDED. INCREASE SIZE'/ &
                5X,'OF MRING. PROGRAM STOPPED IN REPLACE')
  620 FORMAT(///5X,'AT LEAST ONE NEW TETRAHEDRON HAS TOO SMALL A VOLUME' &
               /5X,'PROGRAM STOPPED IN REPLACE')
  640 FORMAT(///5X,'UNABLE TO FIND COMMON FACE BETWEEN ADJACENT CELLS' &
               /5X,'PROGRAM STOPPED IN REPLACE')
  660 FORMAT(//5X,'UNABLE TO FIND EDGE ADDRESS FOR A NEW TETRAHEDRON'/ &
               5X,'PROGRAM STOPPED IN ROUTINE REPLACE') 
  670 FORMAT(//5X,'DIMENSION OF NDG ARRAY EXCEEDED IN ROUTINE REPLACE.'/ &
               5X,'INCREASE SIZE OF MBPTS.')
      END SUBROUTINE REPLACE






!
!     ******************************************************************
!
      SUBROUTINE TRISWP (NPT,N1,N2,N3,J,NFAIL, &
                         X,NDC,NBH,IPROT,NCELL,NDG,IDGP,NDGP,NEDGE, &
                         VOL,XCEN,YCEN,ZCEN,RC,RAT,NPTET,NACPT,TOLV)
!
!     ******************************************************************
!     *                                                                *
!     *  GIVEN TWO TETRAHEDRA WITH COMMON FACE (N1,N2,N3), REMOVE THIS *
!     *  FACE AND INSERT AN EDGE TO CREATE A SET OF THREE TETRAHEDRA.  *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   2001                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: J,NEDGE,NFAIL,NPT,N1,N2,N3
      INTEGER :: NPTET(*),NACPT(*),NDC(4,*),NBH(4,*),IPROT(*),NDG(2,*), &
                 IDGP(*),NDGP(*)
      DOUBLE PRECISION :: TOLV
      DOUBLE PRECISION :: X(3,*),VOL(*),XCEN(*),YCEN(*),ZCEN(*), &
                          RC(*),RAT(*)
 
      INTEGER :: I,IEDG,IM1,IP1,ISMALL,K,KK,KM1,KP1,K1,K2,LA,L1,L2, &
                 M1,M2,M3,M4,NA,NCELL,NOLD
      INTEGER :: IRING(3),NEWC(3),LNBR(3),ISHK(3),MNBR(3),KSHK(3)
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,RADMAX,RATMAX,RCRIN, & 
                          RNX,RNY,RNZ,XCN,YCN,ZCN,VL,VTET1,VTET2,RADC
!
!     ******************************************************************
!
      NFAIL       = 0
!
!     FIND NEIGHBORING TETRAHEDRON LA AND OPPOSITE VERTEX NA.
!
      CALL NEIGHB (J,N1,N2,N3,LA,K1,K2,NDC,NBH)
!
      IF (IPROT(LA).EQ.1) RETURN
      NA          = NDC(1,LA)  +NDC(2,LA)  +NDC(3,LA)  +NDC(4,LA) &
                  -N1  -N2  -N3
      IRING(1)    = N1
      IRING(2)    = N2
      IRING(3)    = N3
!
!     CHECK WHETHER THE TETRAHEDRAL RING FORMS A CONVEX ENSEMBLE
!
      RADMAX  = MAX(RC(J),RC(LA))
      RATMAX  = MAX(RAT(J),RAT(LA))
      DO 10 I=1,3
      IP1     = MOD(I,3)  +1
      IM1     = MOD(I+1,3)  +1
      M1      = IRING(I)
      M2      = IRING(IP1)
      M3      = IRING(IM1)

      CALL CIRCUM (X,M1,M2,NPT,NA,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)

      IF (RADC.GT.RADMAX) GO TO 300
      AREA      = TETAR2 (M1,M2,NPT,NA,X)
      RCRIN     = RADC*AREA/VL
      IF (RCRIN.GT.RATMAX) GO TO 300
      RNX     = COFACT (X(2,M1),X(2,M2),X(2,NA),X(3,M1),X(3,M2),X(3,NA))
      RNY     = COFACT (X(3,M1),X(3,M2),X(3,NA),X(1,M1),X(1,M2),X(1,NA))
      RNZ     = COFACT (X(1,M1),X(1,M2),X(1,NA),X(2,M1),X(2,M2),X(2,NA))
      VTET1   = RNX*(X(1,M3)  -X(1,M1))  +RNY*(X(2,M3)  -X(2,M1)) &
                                         +RNZ*(X(3,M3)  -X(3,M1))
      VTET2   = RNX*(X(1,NPT)  -X(1,M1))  +RNY*(X(2,NPT)  -X(2,M1)) &
                                          +RNZ*(X(3,NPT)  -X(3,M1))
      IF (VTET1*VTET2.LT.0) GO TO 310
      IF (ABS(VTET1).LT.10.*TOLV) GO TO 310
      IF (ABS(VTET2).LT.10.*TOLV) GO TO 310
   10 CONTINUE
!
!     FIND THE SIX NEIGHBORING TETRAHEDRA LNBR(K),MNBR(K), K=1,3
!
      DO 20 K=1,3
      KP1         = MOD(K,3)  +1

      CALL NEIGHB (J,IRING(K),IRING(KP1),NPT,LNBR(K),KK,ISHK(K),NDC,NBH)

      CALL NEIGHB (LA,IRING(K),IRING(KP1),NA,MNBR(K),KK,KSHK(K),NDC,NBH)

   20 CONTINUE
!
!     ASSIGN ADDRESSES OF NEW CELLS
!
      NEWC(1)     = J
      NEWC(2)     = LA
      NCELL       = NCELL  +1
      IF (NCELL.GT.MXCELL) GO TO 360
      NEWC(3)     = NCELL
!
!     CREATE NDC ARRAY FOR NEW CELLS
!
      DO 30 K=1,3
      KP1       = MOD(K,3)  +1
      NDC(1,NEWC(K)) = NPT
      NDC(2,NEWC(K)) = IRING(K)
      NDC(3,NEWC(K)) = IRING(KP1)
      NDC(4,NEWC(K)) = NA
      NACPT(NEWC(K)) = 1
      IPROT(NEWC(K)) = 0
   30 CONTINUE
!
!     COMPUTE VOLUME,CIRCUMCENTER AND CIRCUMRADIUS FOR NEW CELLS
!
      DO 40 K=1,3
      M1        = NDC(1,NEWC(K))
      M2        = NDC(2,NEWC(K))
      M3        = NDC(3,NEWC(K))
      M4        = NDC(4,NEWC(K))

      CALL CIRCUM (X,M1,M2,M3,M4,XCN,YCN,ZCN,VL,RADC,ISMALL,TOLV)


      IF (ISMALL.EQ.1) WRITE (6,910) NPT,NA,(IRING(KK),KK=1,3)
  910 FORMAT('NPT ',I6,' NA ',I6,' IRING ',4I6)
      IF (ISMALL.EQ.1) WRITE (6,960) K,VL, &
                       M1,X(1,M1),X(2,M1),X(3,M1), &
                       M2,X(1,M2),X(2,M2),X(3,M2), &
                       M3,X(1,M3),X(2,M3),X(3,M3), &
                       M4,X(1,M4),X(2,M4),X(3,M4)
  960 FORMAT('RING ',I2,' VOLUME = ',E13.5/ &
             'M1 ',I6,' X = ',F10.4,' Y = ',F10.4,' Z = ',F10.4/ &
             'M2 ',I6,' X = ',F10.4,' Y = ',F10.4,' Z = ',F10.4/ &
             'M3 ',I6,' X = ',F10.4,' Y = ',F10.4,' Z = ',F10.4/ &
             'M4 ',I6,' X = ',F10.4,' Y = ',F10.4,' Z = ',F10.4)

      IF (ISMALL.EQ.1) GO TO 350
      XCEN(NEWC(K)) = XCN
      YCEN(NEWC(K)) = YCN
      ZCEN(NEWC(K)) = ZCN
      VOL(NEWC(K))  = VL
      RC(NEWC(K))   = RADC
      AREA       = TETAR (NEWC(K),X,NDC)
      RAT(NEWC(K)) = RC(NEWC(K))*AREA/VOL(NEWC(K))
   40 CONTINUE
!
!     UPDATE NBH ARRAY
!
      DO 50 K=1,3
      KM1       = MOD(K+1,3)  +1
      KP1       = MOD(K,3)  +1
      NBH(ISHK(K),LNBR(K)) = NEWC(K)
      NBH(KSHK(K),MNBR(K)) = NEWC(K)
      NBH(1,NEWC(K)) = NEWC(KM1)
      NBH(2,NEWC(K)) = NEWC(KP1)
      NBH(3,NEWC(K)) = LNBR(K)
      NBH(4,NEWC(K)) = MNBR(K)
   50 CONTINUE
!
!     INSERT NEW EDGE
!
      L1        = MIN(NPT,NA)
      L2        = MAX(NPT,NA)
      IEDG      = IDGP(L1)
      IF (IEDG.EQ.0) GO TO 60
   55 IF (NDGP(IEDG).EQ.0) GO TO 60
      IEDG      = NDGP(IEDG)
      GO TO 55
   60 NEDGE     = NEDGE  +1
      IF (NEDGE.GT.MXEDGE) GO TO 330
      NDG(1,NEDGE) = L1
      NDG(2,NEDGE) = L2
      NDGP(NEDGE)  = 0
      IF (IEDG.EQ.0) IDGP(L1) = NEDGE
      IF (IEDG.GT.0) NDGP(IEDG) = NEDGE
!
!     UPDATE NPTET ARRAY
!
      NPTET(NPT) = NEWC(1)
      NPTET(NA)  = NEWC(1)
      DO 80 K=1,3
      NPTET(IRING(K)) = NEWC(K)
   80 CONTINUE
      RETURN
  300 NFAIL      = 1
!     WRITE (6,600)
      RETURN
  310 NFAIL      = 1
!     WRITE (6,610)
      RETURN
  320 WRITE (6,620) NOLD,L1,L2,NPT,NA
      STOP
  330 WRITE (6,630)
      STOP
  350 WRITE (6,650)
      STOP
  360 WRITE (6,660)
      STOP
  600 FORMAT(5X,'A NEW TETRAHEDRON IS WORSE THAN AN EXISTING ONE')
  610 FORMAT(5X,'TETRAHEDRAL RING IS NON-CONVEX')
  620 FORMAT(///5X,'EDGE TO BE INSERTED ALREADY EXISTS'/ &
                5X,'EDGE  ',I6,' VERTICES ',2I6,' NPT ',I6,' NA ',I6/ &
              //5X,'PROGRAM STOPPED IN TRISWP')
  630 FORMAT(///5X,'DIMENSION OF ARRAY NDG EXCEEDED IN ROUTINE TRISWP.'/ &
              //5X,'INCREASE SIZE OF MNODE.'/ &
              //5X,'PROGRAM STOPPED IN TRISWP')
  650 FORMAT(///5X,'AT LEAST ONE OF THE NEW TETRAHEDRA HAS TOO'/ &
                5X,'SMALL A VOLUME.'/ &
              //5X,'PROGRAM STOPPED IN TRISWP')
  660 FORMAT(///5X,'DIMENSION OF ARRAY NDC EXCEEDS ALLOWED VALUE.'/ &
              //5X,'INCREASE SIZE OF MCELL.'/ &
              //5X,'PROGRAM STOPPED IN TRISWP')
      END SUBROUTINE TRISWP






!
!     ******************************************************************
!
      SUBROUTINE TESSEL (IRING,NRING,NFAD,NBHKP,IEDKP,NEDK,IMEET,ITOUCH, &
                         X,IPOINT,NPOINT,DX,DY,DZ,DS,NP)
!
!     ******************************************************************
!     *                                                                *
!     *  CREATE FACET TRIANGULATION OF THE DOMAIN INSIDE THE RING      *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IMEET,ITOUCH,NEDK,NRING
      INTEGER :: NPOINT(*),IPOINT(*)
      INTEGER :: IRING(*),NFAD(3,*),NBHKP(3,*),IEDKP(4,*),NP(*)
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: DX(*),DY(*),DZ(*),DS(*)
      
      INTEGER :: I,IA,IB,IMIN,IMIN1,IM1,IP1,K,KPRE,L,LA,LB,NCNT,NEDTOT, &
                 NSWIT,NUM
!
!     ******************************************************************
!
!     LOAD RING EDGES INTO EDGE DATA STRUCTURE
!
      NEDTOT       = 2*NRING  -3
      DO 2 I=1,NEDTOT
      NPOINT(I)    = 0
   2  CONTINUE
      DO 4 I=1,NRING
      IPOINT(IRING(I)) = 0
      NP(I)        = IRING(I)
   4  CONTINUE
      DO 10 I=1,NRING
      IP1          = MOD(I,NRING)  +1
      IA           = MIN(IRING(I),IRING(IP1))
      IB           = MAX(IRING(I),IRING(IP1))
      IEDKP(1,I)   = IA
      IEDKP(2,I)   = IB
      IEDKP(3,I)   = I
      IEDKP(4,I)   = 0
      K            = IPOINT(IA)
   6  IF (K.EQ.0) GO TO 8
      KPRE         = K
      K            = NPOINT(KPRE)
      GO TO 6
   8  IF (IPOINT(IA).NE.0) NPOINT(KPRE) = I
      IF (IPOINT(IA).EQ.0) IPOINT(IA) = I
      NPOINT(I)    = 0
   10 CONTINUE
      NEDK         = NRING
      NUM          = NRING
!
!     ITERATE TO CREATE TRIANGULATION OF RING INTERIOR
!
      NSWIT        = 0
   20 IF (NUM.EQ.3) GO TO 60

      CALL ANGFND (NUM,IRING,IMIN,IMEET,ITOUCH,X,DX,DY,DZ,DS)

      IM1          = MOD(NUM+IMIN-2,NUM)  +1
      IP1          = MOD(IMIN,NUM)  +1
      NSWIT        = NSWIT  +1
      NFAD(1,NSWIT) = IRING(IMIN)
      NFAD(2,NSWIT) = IRING(IP1)
      NFAD(3,NSWIT) = IRING(IM1)
      NEDK         = NEDK  +1
      IA           = MIN(IRING(IM1),IRING(IP1))
      IB           = MAX(IRING(IM1),IRING(IP1))
      IEDKP(1,NEDK) = IA
      IEDKP(2,NEDK) = IB
      IEDKP(3,NEDK) = -NSWIT
      IEDKP(4,NEDK) = 0
      K            = IPOINT(IA)
   22 IF (K.EQ.0) GO TO 24
      KPRE         = K
      K            = NPOINT(KPRE)
      GO TO 22
   24 IF (IPOINT(IA).NE.0) NPOINT(KPRE) = NEDK
      IF (IPOINT(IA).EQ.0) IPOINT(IA) = NEDK
      IA           = MIN(IRING(IM1),IRING(IMIN))
      IB           = MAX(IRING(IM1),IRING(IMIN))
      K            = IPOINT(IA)
   25 IF (IB.EQ.IEDKP(2,K)) GO TO 30
      K            = NPOINT(K)
      IF (K.EQ.0) GO TO 330
      GO TO 25
   30 IEDKP(4,K)   = -NSWIT
      IA           = MIN(IRING(IP1),IRING(IMIN))
      IB           = MAX(IRING(IP1),IRING(IMIN))
      K            = IPOINT(IA)
   35 IF (IB.EQ.IEDKP(2,K)) GO TO 40
      K            = NPOINT(K)
      IF (K.EQ.0) GO TO 340
      GO TO 35
   40 IEDKP(4,K)   = -NSWIT
      NCNT         = 0
      DO 50 I=1,NUM
      IF (I.EQ.IMIN) GO TO 50
      NCNT         = NCNT  +1
      IRING(NCNT)  = IRING(I)
   50 CONTINUE
      NUM          = NCNT
      GO TO 20
!
!     ONLY ONE TRIANGLE REMAINING
!
   60 NSWIT         = NSWIT  +1
      IF (NSWIT.NE.NRING-2) GO TO 320
      NFAD(1,NSWIT) = IRING(1)
      NFAD(2,NSWIT) = IRING(2)
      NFAD(3,NSWIT) = IRING(3)
      IA           = MIN(IRING(1),IRING(2))
      IB           = MAX(IRING(1),IRING(2))
      K            = IPOINT(IA)
   65 IF (IB.EQ.IEDKP(2,K)) GO TO 70
      K            = NPOINT(K)
      GO TO 65
   70 IEDKP(4,K)   = -NSWIT
      IA           = MIN(IRING(2),IRING(3))
      IB           = MAX(IRING(2),IRING(3))
      K            = IPOINT(IA)
   75 IF (IB.EQ.IEDKP(2,K)) GO TO 80
      K            = NPOINT(K)
      GO TO 75
   80 IEDKP(4,K)   = -NSWIT
      IA           = MIN(IRING(3),IRING(1))
      IB           = MAX(IRING(3),IRING(1))
      K            = IPOINT(IA)
   85 IF (IB.EQ.IEDKP(2,K)) GO TO 90
      K            = NPOINT(K)
      GO TO 85
   90 IEDKP(4,K)   = -NSWIT
!
!     ASSEMBLE ADJACENCY INFORMATION FOR RING TRIANGLES
! 
      DO 100 I=1,NRING
      IPOINT(NP(I)) = 0
      IRING(I)     = NP(I)
  100 CONTINUE
      DO 110 K=1,NEDK
      NPOINT(K)    = 0
  110 CONTINUE
      DO 120 K=1,NEDK
      LA           = IABS(IEDKP(3,K))
      LB           = IABS(IEDKP(4,K))
      IF (LA.EQ.0.OR.LB.EQ.0) GO TO 300
      NPOINT(LB)   = NPOINT(LB)  +1
      NBHKP(NPOINT(LB),LB) = IEDKP(3,K)
      IF (IEDKP(3,K).GT.0) GO TO 120
      NPOINT(LA)   = NPOINT(LA)  +1
      NBHKP(NPOINT(LA),LA) = IEDKP(4,K)
  120 CONTINUE
      DO 125 L=1,NSWIT
      IF (NPOINT(L).NE.3) GO TO 310
      NPOINT(L)    = 0
  125 CONTINUE
      RETURN
  300 WRITE (6,600) LA,LB
      STOP
  310 WRITE (6,610)
      STOP
  320 WRITE (6,620) NRING,NSWIT
      STOP
  330 WRITE (6,630) IA,IB
      STOP
  340 WRITE (6,640) IA,IB
      STOP
  600 FORMAT(/5X,'A ZERO ADDRESS HAS BEEN FOUND FOR A TRIANGLE', &
              5X,'ADDRESSES ARE ',I6,' AND ',I6, &
              5X,'PROGRAM STOPPED IN ROUTINE TESSEL')
  610 FORMAT(/5X,'INCORRECT ASSEMBLY OF TRIANGLES', &
              5X,'PROGRAM STOPPED IN ROUTINE TESSEL')
  620 FORMAT(/5X,'NRING IS NOT EQUAL TO NSWIT PLUS TWO', &
              5X,'NRING = ',I3,'  NSWIT = ',I3, &
              5X,'PROGRAM STOPPED IN ROUTINE TESSEL')
  630 FORMAT(/5X,'EDGE WITH ZERO ADDRESS FOUND AFTER LABEL 25', &
              5X,'IA = ',I6,' IB = ',I6, &
              5X,'PROGRAM STOPPED IN ROUTINE TESSEL')
  640 FORMAT(/5X,'EDGE WITH ZERO ADDRESS FOUND AFTER LABEL 35', &
              5X,'IA = ',I6,' IB = ',I6, &
              5X,'PROGRAM STOPPED IN ROUTINE TESSEL')
      END SUBROUTINE TESSEL








!
!     ******************************************************************
!
      SUBROUTINE ANGFND (NUM,IRING,IMIN,IMEET,ITOUCH,X,DX,DY,DZ,DS)
!
!     ******************************************************************
!     *                                                                *
!     *  SEARCH THROUGH RING OF ADJACENT EDGES TO FIND VERTEX WHICH    *
!     *  HAS THE SMALLEST INTERIOR ANGLE                               *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IMEET,IMIN,ITOUCH,NUM
      INTEGER :: IRING(*)
      DOUBLE PRECISION :: X(3,*)
      DOUBLE PRECISION :: DX(*),DY(*),DZ(*),DS(*)
      
      INTEGER :: I,IM1,IP1
      DOUBLE PRECISION :: ANG,ANGMIN,BETA,CANG,DET,DIV,DS1SQ,DS2SQ,FAC, &
                          FACT1,FACT2,GAMMA,PI,PROD,RNX,RPX,RNY,RPY, &
                          RNZ,RPZ,SANG,SUM
!
!     ******************************************************************
!
      PI         = 4*ATAN(1.)
      SUM        = 0.
      IMIN       = 0
      DO 10 I=1,NUM
      IM1        = MOD(NUM+I-2,NUM)  +1
      DX(I)      = X(1,IRING(IM1))  -X(1,IRING(I))
      DY(I)      = X(2,IRING(IM1))  -X(2,IRING(I))
      DZ(I)      = X(3,IRING(IM1))  -X(3,IRING(I))
      DS(I)      = 1./SQRT(DX(I)*DX(I)  +DY(I)*DY(I)  +DZ(I)*DZ(I))
   10 CONTINUE
      RPX        = 0.
      RPY        = 0.
      RPZ        = 0.
      DO 20 I=1,NUM
      IP1        = MOD(I,NUM)  +1
      FAC        = DS(I)*DS(IP1)
      RPX        = RPX  -FAC*(DY(IP1)*DZ(I)  -DY(I)*DZ(IP1))
      RPY        = RPY  -FAC*(DZ(IP1)*DX(I)  -DZ(I)*DX(IP1)) 
      RPZ        = RPZ  -FAC*(DX(IP1)*DY(I)  -DX(I)*DY(IP1))
   20 CONTINUE
      FAC        = 1./SQRT(RPX*RPX  +RPY*RPY  +RPZ*RPZ)
      RPX        = FAC*RPX
      RPY        = FAC*RPY
      RPZ        = FAC*RPZ
      DO 30 I=1,NUM
      IP1        = MOD(I,NUM)  +1
      IM1        = MOD(NUM+I-2,NUM)  +1
      PROD       = -DX(I)*DX(IP1)  -DY(I)*DY(IP1)  -DZ(I)*DZ(IP1)
      FACT1      = DS(I)*DS(IP1)  -PROD
      FACT2      = DS(I)*DS(IP1)  +PROD
      DET        = 1./(FACT1*FACT2)
      DS1SQ      = DS(I)*DS(I)
      DS2SQ      = DS(IP1)*DS(IP1)
      BETA       = DS2SQ*(RPX*DX(I)  +RPY*DY(I)  +RPZ*DZ(I)) &
                  +PROD*(RPX*DX(IP1)  +RPY*DY(IP1)  +RPZ*DZ(IP1))*DET
      GAMMA      = DS1SQ*(RPX*DX(IP1)  +RPY*DY(IP1)  +RPZ*DZ(IP1)) &
                  +PROD*(RPX*DX(I)  +RPY*DY(I)  +RPZ*DZ(I))*DET
      RNX        = RPX  -BETA*DX(I)  -GAMMA*DX(IP1)
      RNY        = RPY  -BETA*DY(I)  -GAMMA*DY(IP1)
      RNZ        = RPZ  -BETA*DZ(I)  -GAMMA*DZ(IP1)
      FAC        = 1./SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      RNX        = FAC*RNX
      RNY        = FAC*RNY
      RNZ        = FAC*RNZ
      DIV        = DS(I)*DS(IP1)
      CANG       = DIV*PROD
      SANG       = -DIV*((DY(IP1)*DZ(I)  -DY(I)*DZ(IP1))*RNX &
                        +(DZ(IP1)*DX(I)  -DZ(I)*DX(IP1))*RNY &
                        +(DX(IP1)*DY(I)  -DX(I)*DY(IP1))*RNZ)
      ANG        = ATAN2(SANG,CANG)
      IF (ANG.LT.0.) ANG = ANG  +2.*PI
      IF (IMIN.EQ.0) GO TO 25
      IF (ANG.GT.ANGMIN) GO TO 30
      IF (IRING(I).EQ.IMEET.AND.IRING(IM1).EQ.ITOUCH) GO TO 25
      IF (IRING(I).EQ.IMEET.AND.IRING(IP1).EQ.ITOUCH) GO TO 25
      IF (IRING(I).EQ.ITOUCH.AND.IRING(IM1).EQ.IMEET) GO TO 25
      IF (IRING(I).EQ.ITOUCH.AND.IRING(IP1).EQ.IMEET) GO TO 25
      IF (IRING(I).EQ.IMEET.OR.IRING(I).EQ.ITOUCH) GO TO 30
   25 ANGMIN     = ANG
      IMIN       = I
   30 CONTINUE
      RETURN
      END SUBROUTINE ANGFND






!
!     ******************************************************************
!
      SUBROUTINE SLIVER (IRING,NRING,NA,NB,NFAIL,VORIG, &
                         X,DENS,ITYP,IPOINT,NNODE, &
                         NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR,VLT)
!
!     ******************************************************************
!     *                                                                *
!     *  INSERT A POINT INSIDE TETRAHEDRAL RING IN ORDER TO REMOVE     *
!     *  SLIVER                                                        *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IOCTR,NA,NB,NFAIL,NNODE,NRING
      INTEGER :: IPOINT(*),IRING(*),ITYP(*),NREF(*),NLINK(*),NOCTR(2,*)
      DOUBLE PRECISION :: VORIG
      DOUBLE PRECISION :: X(3,*),DENS(*)
      DOUBLE PRECISION :: XFAR(2),YFAR(2),ZFAR(2),VLT(*)

      INTEGER :: I,IM1,M1,M2,N,NCNT,NCYC,NRM1
      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,DENSUM,FAC,RNX,RNY,RNZ,TOL, &
                          VTET1,VTET2,XPT,YPT,ZPT
!
!     ******************************************************************
!
      TOL       = 1.E-15
      XPT       = 0.
      YPT       = 0.
      ZPT       = 0.
      DENSUM    = 0.
      NRM1      = NRING  -1
      DO 10 I=2,NRM1
      XPT       = XPT  +X(1,IRING(I))
      YPT       = YPT  +X(2,IRING(I))
      ZPT       = ZPT  +X(3,IRING(I))
      DENSUM    = DENSUM  +DENS(IRING(I))
   10 CONTINUE
      FAC       = 1./FLOAT(NRING-2)
      XPT       = .5*FAC*XPT  +.25*(X(1,NA)  +X(1,NB))
      YPT       = .5*FAC*YPT  +.25*(X(2,NA)  +X(2,NB))
      ZPT       = .5*FAC*ZPT  +.25*(X(3,NA)  +X(3,NB))
      DENSUM    = .5*FAC*DENSUM  +.25*(DENS(NA)  +DENS(NB))
      NCYC      = 0
   12 NCNT      = 0
      DO 20 I=1,NRING
      IM1       = MOD(NRING-2+I,NRING)  +1
      M1        = IRING(I)
      M2        = IRING(IM1)
      RNX     = COFACT (X(2,M1),X(2,M2),X(2,NA),X(3,M1),X(3,M2),X(3,NA))
      RNY     = COFACT (X(3,M1),X(3,M2),X(3,NA),X(1,M1),X(1,M2),X(1,NA))
      RNZ     = COFACT (X(1,M1),X(1,M2),X(1,NA),X(2,M1),X(2,M2),X(2,NA))
      VTET2     = RNX*(XPT  -X(1,M1))  +RNY*(YPT  -X(2,M1)) &
                                       +RNZ*(ZPT  -X(3,M1))
      NCNT      = NCNT  +1
      VLT(NCNT) = VTET2
      IF (ABS(VTET2).LT.VORIG) GO TO 210
      IF (I.EQ.1) GO TO 15
      IF (VTET1*VTET2.LE.TOL) GO TO 100
      GO TO 20
   15 VTET1     = VTET2
   20 CONTINUE
      DO 30 I=1,NRING
      IM1       = MOD(NRING-2+I,NRING)  +1
      M1        = IRING(I)
      M2        = IRING(IM1)
      RNX     = COFACT (X(2,M1),X(2,M2),X(2,NB),X(3,M1),X(3,M2),X(3,NB))
      RNY     = COFACT (X(3,M1),X(3,M2),X(3,NB),X(1,M1),X(1,M2),X(1,NB))
      RNZ     = COFACT (X(1,M1),X(1,M2),X(1,NB),X(2,M1),X(2,M2),X(2,NB))
      VTET2     = RNX*(XPT  -X(1,M1))  +RNY*(YPT  -X(2,M1)) &
                                       +RNZ*(ZPT  -X(3,M1))
      NCNT      = NCNT  +1
      VLT(NCNT) = VTET2
      IF (ABS(VTET2).LT.VORIG) GO TO 210
      IF (VTET1*VTET2.GE.TOL) GO TO 100
   30 CONTINUE
      NNODE     = NNODE  +1
      IF (NNODE.GT.MXNODE) GO TO 230
      N         = NNODE
      X(1,N)    = XPT
      X(2,N)    = YPT
      X(3,N)    = ZPT
      DENS(N)   = DENSUM
      ITYP(N)   = 8
      IPOINT(N) = 0

      CALL OCTFIL (N,X,NOCTR,IOCTR,NLINK,NREF,XFAR,YFAR,ZFAR)

      RETURN
  100 IF (NCYC.EQ.2) GO TO 200
      NCYC      = NCYC  +1
      XPT       = .5*XPT  +.25*(X(1,NA)  +X(1,NB))
      YPT       = .5*YPT  +.25*(X(2,NA)  +X(2,NB))
      ZPT       = .5*ZPT  +.25*(X(3,NA)  +X(3,NB))
      DENSUM    = .5*DENSUM  +.25*(DENS(NA)  +DENS(NB))
      GO TO 12
  200 NFAIL     = 2
!     WRITE (6,600)
      RETURN
  210 NFAIL     = 2
!     WRITE (6,610) VORIG,VTET2
      RETURN
  230 WRITE (6,630)
      STOP
  600 FORMAT(//'POINT TO BE INSERTED LIES OUTSIDE TETRAHEDRAL RING')
  610 FORMAT(//'A NEW TETRAHEDRON IS TOO SMALL'/ &
               'VORIG = ',E13.5,'   VTET = ',E13.5) 
  630 FORMAT(//5X,'NNODE EXCEEDS DIMENSION OF ARRAY X.'/ &
               5X,'INCREASE SIZE OF MNODE.'/ &
               5X,'PROGRAM STOPPED IN ROUTINE SLIVER')
      END SUBROUTINE SLIVER






!
!     ******************************************************************
!
      SUBROUTINE CIRCUM (X,N1,N2,N3,N4,XCN,YCN,ZCN,VTET,RCIR,IFLAG,TOL)
!
!     ******************************************************************
!     *                                                                *
!     *  FIND VOLUME, CIRCUMRADIUS AND CIRCUMCENTER COORDINATES OF     *
!     *  TETRAHEDRON (N1,N2,N3,N4).                                    *
!     *  IFLAG = 0 FOR IF VOLUME EXCEEDS THRESHOLD.                    *
!     *  IFLAG = 1 IF VOLUME IS TOO SMALL.                             *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: IFLAG,N1,N2,N3,N4
      DOUBLE PRECISION :: XCN,YCN,ZCN,VTET,RCIR,TOL         
      DOUBLE PRECISION :: X(3,*)
      
      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,C5,C6,C9,DET,FAC,H1,H2,H3, &
                          RNX,RNY,RNZ,VCELL,XSHF,YSHF,ZSHF   
!
!     ******************************************************************
!
      IFLAG     = 0
!
!     COMPUTE VOLUME OF TETRAHEDRON
!
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      DET     = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                         +RNZ*(X(3,N4)  -X(3,N1))
      VTET    = ABS(DET)
!
!     COMPUTE CIRCUMCENTER COORDINATES
!
      RCIR      = 1.E10
      IF (VTET.LT.TOL) GO TO 30
      FAC       = .5/DET
      XSHF      = .5*(X(1,N1)  +X(1,N2)  +X(1,N3)  +X(1,N4))
      YSHF      = .5*(X(2,N1)  +X(2,N2)  +X(2,N3)  +X(2,N4))
      ZSHF      = .5*(X(3,N1)  +X(3,N2)  +X(3,N3)  +X(3,N4))
      H1        = (X(1,N4)  -X(1,N1))*(X(1,N4)  -XSHF  +X(1,N1)) &
                 +(X(2,N4)  -X(2,N1))*(X(2,N4)  -YSHF  +X(2,N1)) &
                 +(X(3,N4)  -X(3,N1))*(X(3,N4)  -ZSHF  +X(3,N1))
      H2        = (X(1,N3)  -X(1,N1))*(X(1,N3)  -XSHF  +X(1,N1)) &
                 +(X(2,N3)  -X(2,N1))*(X(2,N3)  -YSHF  +X(2,N1)) &
                 +(X(3,N3)  -X(3,N1))*(X(3,N3)  -ZSHF  +X(3,N1))
      H3        = (X(1,N2)  -X(1,N1))*(X(1,N2)  -XSHF  +X(1,N1)) &
                 +(X(2,N2)  -X(2,N1))*(X(2,N2)  -YSHF  +X(2,N1)) &
                 +(X(3,N2)  -X(3,N1))*(X(3,N2)  -ZSHF  +X(3,N1))
      C5        = H2*(X(3,N2)  -X(3,N1))  -H3*(X(3,N3)  -X(3,N1))
      C6        = H2*(X(2,N2)  -X(2,N1))  -H3*(X(2,N3)  -X(2,N1))
      C9        = H3*(X(1,N3)  -X(1,N1))  -H2*(X(1,N2)  -X(1,N1))
      XCN       = .5*XSHF  +(H1*RNX  +(X(2,N4)  -X(2,N1))*C5 &
                                     -(X(3,N4)  -X(3,N1))*C6)*FAC
      YCN       = .5*YSHF  +(-(X(1,N4)  -X(1,N1))*C5  +H1*RNY &
                             -(X(3,N4)  -X(3,N1))*C9)*FAC 
      ZCN       = .5*ZSHF  +((X(1,N4)  -X(1,N1))*C6 &
                            +(X(2,N4)  -X(2,N1))*C9  +H1*RNZ)*FAC
!
!     COMPUTE CIRCUMRADIUS
!
      RCIR      = SQRT((X(1,N1)  -XCN)**2  +(X(2,N1)  -YCN)**2 &
                                           +(X(3,N1)  -ZCN)**2)
      RETURN
   30 VCELL     = VTET/6.
!     WRITE (6,600) VCELL
      IFLAG     = 1
      RETURN
  600 FORMAT(//5X,'TETRAHEDRON WITH AN EXTREMELY SMALL VOLUME FOUND'// &
                  ' IN ROUTINE CIRCUM'//5X,'VOLUME = ',E13.5/)
      END SUBROUTINE CIRCUM






!
!     ******************************************************************
!
      SUBROUTINE CIRCUM2 (X,N1,N2,N3,XPT,YPT,ZPT, &
                          XCN,YCN,ZCN,VTET,RCIR,IFLAG,TOL)
!
!     ******************************************************************
!     *                                                                *
!     *  FIND VOLUME, CIRCUMRADIUS AND CIRCUMCENTER COORDINATES OF     *
!     *  TETRAHEDRON WHOSE BASE HAS THE VERTEX ADDRESSES N1,N2,N3 AND  *
!     *  WHOSE FOURTH VERTEX IS THE POINT (XPT,YPT,ZPT).               *
!     *  IFLAG = 0 FOR IF VOLUME EXCEEDS THRESHOLD.                    *
!     *  IFLAG = 1 IF VOLUME IS TOO SMALL.                             *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: N1,N2,N3
      DOUBLE PRECISION :: XCN,XPT,YCN,YPT,ZCN,ZPT,VTET,RCIR,TOL
      DOUBLE PRECISION :: X(3,*)      

      INTEGER :: IFLAG
      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,C5,C6,C9,DET,FAC,H1,H2,H3, & 
                          RNX,RNY,RNZ,VCELL,XSHF,YSHF,ZSHF
!
!     ******************************************************************
!
      IFLAG     = 0
!
!     COMPUTE VOLUME OF TETRAHEDRON
!
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      DET     = RNX*(XPT  -X(1,N1))  +RNY*(YPT  -X(2,N1)) &
                                     +RNZ*(ZPT  -X(3,N1))
      VTET    = ABS(DET)
!
!     COMPUTE CIRCUMCENTER COORDINATES
!
      RCIR      = 1.E10
      IF (VTET.LT.TOL) GO TO 30
      FAC       = .5/DET
      XSHF      = .5*(X(1,N1)  +X(1,N2)  +X(1,N3)  +XPT)
      YSHF      = .5*(X(2,N1)  +X(2,N2)  +X(2,N3)  +YPT)
      ZSHF      = .5*(X(3,N1)  +X(3,N2)  +X(3,N3)  +ZPT)
      H1        = (XPT  -X(1,N1))*(XPT  -XSHF  +X(1,N1)) &
                 +(YPT  -X(2,N1))*(YPT  -YSHF  +X(2,N1)) &
                 +(ZPT  -X(3,N1))*(ZPT  -ZSHF  +X(3,N1))
      H2        = (X(1,N3)  -X(1,N1))*(X(1,N3)  -XSHF  +X(1,N1)) &
                 +(X(2,N3)  -X(2,N1))*(X(2,N3)  -YSHF  +X(2,N1)) &
                 +(X(3,N3)  -X(3,N1))*(X(3,N3)  -ZSHF  +X(3,N1)) 
      H3        = (X(1,N2)  -X(1,N1))*(X(1,N2)  -XSHF  +X(1,N1)) &
                 +(X(2,N2)  -X(2,N1))*(X(2,N2)  -YSHF  +X(2,N1)) &
                 +(X(3,N2)  -X(3,N1))*(X(3,N2)  -ZSHF  +X(3,N1))
      C5        = H2*(X(3,N2)  -X(3,N1))  -H3*(X(3,N3)  -X(3,N1))
      C6        = H2*(X(2,N2)  -X(2,N1))  -H3*(X(2,N3)  -X(2,N1))
      C9        = H3*(X(1,N3)  -X(1,N1))  -H2*(X(1,N2)  -X(1,N1))
      XCN       = .5*XSHF  +(H1*RNX  +(YPT  -X(2,N1))*C5 &
                                     -(ZPT  -X(3,N1))*C6)*FAC
      YCN       = .5*YSHF  +(-(XPT  -X(1,N1))*C5  +H1*RNY &
                             -(ZPT  -X(3,N1))*C9)*FAC
      ZCN       = .5*ZSHF  +((XPT  -X(1,N1))*C6 &
                            +(YPT  -X(2,N1))*C9  +H1*RNZ)*FAC
!
!     COMPUTE CIRCUMRADIUS
!
      RCIR      = SQRT((XPT  -XCN)**2  +(YPT  -YCN)**2 &
                                       +(ZPT  -ZCN)**2)
      RETURN
   30 VCELL     = VTET/6.
!     WRITE (6,600) VCELL
      IFLAG     = 1
      RETURN
  600 FORMAT(//5X,'TETRAHEDRON WITH AN EXTREMELY SMALL VOLUME FOUND'// &
                  ' IN ROUTINE CIRCUM2'//5X,'VOLUME = ',E13.5/)
      END SUBROUTINE CIRCUM2







!
!     ******************************************************************
!
      SUBROUTINE NEIGHB (L1,N1,N2,N3,L2,K1,K2,NDC,NBH)
!
!     ******************************************************************
!     *                                                                *
!     *  GIVEN FACE (N1,N2,N3) OF TETRAHEDRON L1, FIND THE NEIGHBORING *
!     *  TETRAHEDRON, L2, WHICH SHARES THIS FACE.                      *
!     *  LABELS K1 AND K2 ARE INDICES OF ARRAY NBH.                    *
!     *  THUS, NBH(K1,L1) = L2 AND NBH(K2,L2) = L1                     *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: K1,K2,L1,L2,N1,N2,N3
      INTEGER :: NDC(4,*),NBH(4,*)

      INTEGER :: K,KK,M,MEND,MMAX,MMIN,MSUM,M1,M2,M3,M4,NMAX,NMIN,NSUM
!
!     ******************************************************************
!
      NMIN      = MIN(N1,N2,N3)
      NMAX      = MAX(N1,N2,N3)
      NSUM      = N1  +N2  +N3
      DO 20 K=1,4
      K1        = K
      L2        = NBH(K1,L1)
      M1        = NDC(1,L2)
      M2        = NDC(2,L2)
      M3        = NDC(3,L2)
      M4        = NDC(4,L2)
      MEND      = M4
   10 MMIN      = MIN(M1,M2,M3)
      MMAX      = MAX(M1,M2,M3)
      MSUM      = M1  +M2  +M3
      IF (MMIN.EQ.NMIN.AND.MMAX.EQ.NMAX.AND.MSUM.EQ.NSUM) GO TO 30
      IF (M1.EQ.MEND) GO TO 20
      M         = M1
      M1        = M2
      M2        = M3
      M3        = M4
      M4        = M
      GO TO 10
   20 CONTINUE
      GO TO 110
!
!     NEIGHBORING TETRAHEDRON L2 HAS BEEN FOUND
!
   30 DO 40 K=1,4
      K2        = K
      IF (NBH(K2,L2).EQ.L1) RETURN
   40 CONTINUE
      WRITE (6,600)
      STOP
  110 WRITE (6,610)

      WRITE (6,900) N1,N2,N3, &
                    L1,(NDC(KK,L1),KK=1,4),L2,(NDC(KK,L2),KK=1,4)
  900 FORMAT('N1,N2,N3 ',3I6,' L1 ',I6,' VERTS ',4I6/ &
                         27X,' L2 ',I6,' VERTS ',4I6)

      STOP
  600 FORMAT(//5X,'UNABLE TO FIND NBH INDEX OF NEIGHBORING TETRAHEDRON', &
               5X,'TO ORIGINAL TETRAHEDRON. PROGRAM STOPPED IN NEIGHB')
  610 FORMAT(//5X,'UNABLE TO FIND A NEIGHBORING TETRAHEDRON THAT HAS', &
               5X,'THE REQUIRED FACE. PROGRAM STOPPED IN NEIGHB')
      END SUBROUTINE NEIGHB








!
!     ******************************************************************
!
      SUBROUTINE LOCK (L1,L2,N1,N2,N3,N4,N5,K1,K2,NDC,NBH)
!
!     ******************************************************************
!     *                                                                *
!     *  FIND COMMON FACE BETWEEN TETRAHEDRA L1 AND L2. THEN SET THE   *
!     *  LABELS N1,N2,N3,N4 TO BE THE VERTEX ADDRESSES OF L1 AND SET   *
!     *  N1,N2,N3,N5 TO BE THE VERTEX ADDRESSES OF L2. LABELS K1 AND   *
!     *  K2 ARE INDICES OF ARRAY NBH.                                  *
!     *  THUS, NBH(K1,L1) = L2 AND NBH(K2,L2) = L1                     *
!     *                                                                *
!     ******************************************************************
!    ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: K1,K2,L1,L2,N1,N2,N3,N4,N5
      INTEGER :: NDC(4,*),NBH(4,*)
      
      INTEGER :: M,MEND,MMAX,MMIN,MSUM,M1,M2,M3,M4,N,NEND,NMAX,NMIN,NSUM
!
!     ******************************************************************
!
      N1        = NDC(1,L1)
      N2        = NDC(2,L1)
      N3        = NDC(3,L1)
      N4        = NDC(4,L1)
      NEND      = N4
      M1        = NDC(1,L2)
      M2        = NDC(2,L2)
      M3        = NDC(3,L2)
      M4        = NDC(4,L2)
   10 MEND      = M4
      NMIN      = MIN(N1,N2,N3)
      NMAX      = MAX(N1,N2,N3)
      NSUM      = N1  +N2  +N3
   20 MMIN      = MIN(M1,M2,M3)
      MMAX      = MAX(M1,M2,M3)
      MSUM      = M1  +M2  +M3
      IF (MMIN.EQ.NMIN.AND.MMAX.EQ.NMAX.AND.MSUM.EQ.NSUM) GO TO 40
      IF (M1.EQ.MEND) GO TO 30
      M         = M1
      M1        = M2
      M2        = M3
      M3        = M4
      M4        = M
      GO TO 20
   30 IF (N1.EQ.NEND) GO TO 100
      N         = N1
      N1        = N2
      N2        = N3
      N3        = N4
      N4        = N
      GO TO 10
!
!     COMMON FACE (N1,N2,N3) HAS BEEN FOUND
!
   40 N5        = M4
      K1        = 0
      IF (NBH(1,L1).EQ.L2) K1 = 1
      IF (NBH(2,L1).EQ.L2) K1 = 2
      IF (NBH(3,L1).EQ.L2) K1 = 3
      IF (NBH(4,L1).EQ.L2) K1 = 4
      K2        = 0
      IF (NBH(1,L2).EQ.L1) K2 = 1
      IF (NBH(2,L2).EQ.L1) K2 = 2
      IF (NBH(3,L2).EQ.L1) K2 = 3
      IF (NBH(4,L2).EQ.L1) K2 = 4
      RETURN
  100 WRITE (6,600)
      STOP
  600 FORMAT(//5X,'UNABLE TO FIND THE COMMON FACE IN ROUTINE LOCK')
      END SUBROUTINE LOCK







!
!     ******************************************************************
!
      FUNCTION FACEAR (X,N1,N2,N3)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPUTE AREA OF FACE J                                        *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      DOUBLE PRECISION :: FACEAR

      INTEGER :: N1,N2,N3
      DOUBLE PRECISION :: AX,AY,AZ
      DOUBLE PRECISION :: X(3,*)
!
!     ******************************************************************
!
      AX        = (X(2,N2)  -X(2,N1))*(X(3,N3)  -X(3,N1)) &
                 -(X(3,N2)  -X(3,N1))*(X(2,N3)  -X(2,N1))
      AY        = (X(3,N2)  -X(3,N1))*(X(1,N3)  -X(1,N1)) &
                 -(X(1,N2)  -X(1,N1))*(X(3,N3)  -X(3,N1))
      AZ        = (X(1,N2)  -X(1,N1))*(X(2,N3)  -X(2,N1)) &
                 -(X(2,N2)  -X(2,N1))*(X(1,N3)  -X(1,N1))
      FACEAR  = 0.5D0*SQRT(AX*AX  +AY*AY  +AZ*AZ)
      RETURN
      END FUNCTION FACEAR






!
!     ******************************************************************
!
      FUNCTION DIHED (N1,N2,N3,N4,X)
!
!     ******************************************************************
!     *                                                                *
!     *  CALCULATE DIHEDRAL ANGLE BETWEEN PLANES DEFINED BY (N1,N2,N3) *
!     *  AND (N2,N1,N4) WHICH INTERSECT ALONG COMMON EDGE (N1,N2)      *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      DOUBLE PRECISION :: DIHED
      
      INTEGER N1,N2,N3,N4
      DOUBLE PRECISION :: X(3,*)
      
      DOUBLE PRECISION :: ASIZ,AX,AY,AZ,A1,A2,A3, &
                          BSIZ,BX,BY,BZ,B1,B2,B3, &
                          COSANG,FAC,PI,PX,PY,PZ,RAD, &
                          SINANG,SX,SY,SZ    
!
!     ******************************************************************
!
      PI      = 4.0D0*ATAN(1.0D0)
      RAD     = 180.0D0/PI
      DIHED   = 0.D0
!
!     CALCULATE UNIT NORMAL (AX,AY,AZ) TO TRIANGLE (N1,N2,N3)
!
      AX      = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      AY      = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      AZ      = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      ASIZ    = SQRT(AX*AX  +AY*AY  +AZ*AZ)
      IF (ASIZ.LT.1.0D-6) RETURN
      FAC     = 1./ASIZ
      AX      = AX*FAC
      AY      = AY*FAC
      AZ      = AZ*FAC
!
!     CALCULATE UNIT NORMAL (BX,BY,BZ) TO TRIANGLE (N2,N1,N4)
!
      BX      = COFACT (X(2,N1),X(2,N4),X(2,N2),X(3,N1),X(3,N4),X(3,N2))
      BY      = COFACT (X(3,N1),X(3,N4),X(3,N2),X(1,N1),X(1,N4),X(1,N2))
      BZ      = COFACT (X(1,N1),X(1,N4),X(1,N2),X(2,N1),X(2,N4),X(2,N2))
      BSIZ    = SQRT(BX*BX  +BY*BY  +BZ*BZ)
      IF (BSIZ.LT.1.0D-6) RETURN
      FAC     = 1./BSIZ
      BX      = BX*FAC
      BY      = BY*FAC
      BZ      = BZ*FAC
!
!     FORM COSINE OF ANGLE BETWEEN THE NORMALS
!
      COSANG  = AX*BX  +AY*BY  +AZ*BZ
!
!     COMPUTE UNIT VECTOR (PX,PY,PZ) DIRECTED FROM POINT N1 TO POINT N2
!
      PX      = X(1,N2)  -X(1,N1)
      PY      = X(2,N2)  -X(2,N1)
      PZ      = X(3,N2)  -X(3,N1)
      FAC     = 1.0D0/SQRT(PX*PX  +PY*PY  +PZ*PZ)
      PX      = PX*FAC
      PY      = PY*FAC
      PZ      = PZ*FAC
!
!     FORM SCALAR TRIPLE PRODUCT OF VECTORS (A,B,P) TO OBTAIN THE SINE
!     OF THE ANGLE BETWEEN THE NORMALS
!
      SX      = AY*BZ  -AZ*BY
      SY      = AZ*BX  -AX*BZ
      SZ      = AX*BY  -AY*BX
      SINANG  = SX*PX  +SY*PY  +SZ*PZ
      DIHED   = 180.D0 +RAD*ATAN2(SINANG,COSANG)
      IF (DIHED.GE.360.D0) DIHED = DIHED  -360.D0
      RETURN
      END FUNCTION DIHED






!
!     ******************************************************************
!
      FUNCTION TETAR (L,X,NDC)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPUTE COMBINED AREA OF THE FOUR FACES OF TETRAHEDRON L      *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      DOUBLE PRECISION :: TETAR

      INTEGER :: L
      INTEGER :: NDC(4,*)
      DOUBLE PRECISION :: X(3,*)

      INTEGER :: N,NEND,N1,N2,N3,N4
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,RNX,RNY,RNZ
!
!     ******************************************************************
! 
      N1      = NDC(1,L)
      N2      = NDC(2,L)
      N3      = NDC(3,L)
      N4      = NDC(4,L)
      NEND    = N4
      AREA    = 0.
   10 RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      AREA    = AREA  +SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      IF (N1.EQ.NEND) GO TO 20
      N       = N1
      N1      = N2
      N2      = N3
      N3      = N4
      N4      = N
      GO TO 10
   20 TETAR   = AREA
      RETURN
      END FUNCTION TETAR






!
!     ******************************************************************
!
      FUNCTION TETAR2 (M1,M2,M3,M4,X)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPUTE COMBINED AREA OF THE FOUR FACES OF THE TETRAHEDRON    *
!     *  WHOSE VERTEX ADDRESSES ARE N1, N2, N3 AND N4                  *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      DOUBLE PRECISION :: TETAR2

      INTEGER :: M1,M2,M3,M4
      DOUBLE PRECISION :: X(3,*)
      
      INTEGER :: N,NEND,N1,N2,N3,N4
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,RNX,RNY,RNZ
!
!     ******************************************************************
!
      N1      = M1
      N2      = M2
      N3      = M3
      N4      = M4
      NEND    = N4
      AREA    = 0.
   10 RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      AREA    = AREA  +SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      IF (N1.EQ.NEND) GO TO 20
      N       = N1
      N1      = N2
      N2      = N3
      N3      = N4
      N4      = N
      GO TO 10
   20 TETAR2  = AREA
      RETURN
      END FUNCTION TETAR2






!
!     ******************************************************************
!
      FUNCTION TETAR3 (M1,M2,M3,XPT,YPT,ZPT,X)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPUTE COMBINED AREA OF THE FOUR FACES OF THE TETRAHEDRON    *
!     *  WHOSE BASE HAS THE VERTEX ADDRESSES N1, N2, N3 AND WHOSE      *
!     *  FOURTH VERTEX IS GIVEN BY THE POINT (XPT,YPT,ZPT).            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      DOUBLE PRECISION :: TETAR3

      INTEGER :: M1,M2,M3
      DOUBLE PRECISION :: XPT,YPT,ZPT
      DOUBLE PRECISION :: X(3,*)
      
      INTEGER :: N,NEND,N1,N2,N3
      DOUBLE PRECISION :: AREA,A1,A2,A3,B1,B2,B3,RNX,RNY,RNZ
!
!     ******************************************************************
!
      N1      = M1
      N2      = M2
      N3      = M3
      NEND    = N3
      AREA    = 0.
   10 RNX     = COFACT (X(2,N1),X(2,N2),YPT,X(3,N1),X(3,N2),ZPT)
      RNY     = COFACT (X(3,N1),X(3,N2),ZPT,X(1,N1),X(1,N2),XPT)
      RNZ     = COFACT (X(1,N1),X(1,N2),XPT,X(2,N1),X(2,N2),YPT)
      AREA    = AREA  +SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      IF (N1.EQ.NEND) GO TO 20
      N       = N1
      N1      = N2
      N2      = N3
      N3      = N
      GO TO 10
   20 RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      TETAR3  = AREA  +SQRT(RNX*RNX  +RNY*RNY  +RNZ*RNZ)
      RETURN
      END FUNCTION TETAR3






!
!     ******************************************************************
!
      SUBROUTINE FANGLE (J,X,NFCE,ANGL1,ANGL2,ANGL3,Q)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPUTE FACE ANGLES AND QUALITY MEASURE FOR FACE J            *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: J
      INTEGER :: NFCE(3,*)
      DOUBLE PRECISION :: ANGL1,ANGL2,ANGL3,Q
      DOUBLE PRECISION :: X(3,*)
      
      INTEGER ::N1,N2,N3
      DOUBLE PRECISION :: CHALF,C1,C2,C3,DA,DB,DC,PROD1,PROD2,PROD3, & 
                          RAD,S1,S2,S3,T1,T2,T3
!
!     ******************************************************************
!
      RAD       = 45.0D0/ATAN(1.0D0)
      N1        = NFCE(1,J)
      N2        = NFCE(2,J)
      N3        = NFCE(3,J)
      PROD1     = (X(1,N2)  -X(1,N1))*(X(1,N3)  -X(1,N1)) &
                 +(X(2,N2)  -X(2,N1))*(X(2,N3)  -X(2,N1)) &
                 +(X(3,N2)  -X(3,N1))*(X(3,N3)  -X(3,N1))
      DA        = SQRT((X(1,N2)  -X(1,N1))**2  +(X(2,N2)  -X(2,N1))**2 &
                                               +(X(3,N2)  -X(3,N1))**2)
      DB        = SQRT((X(1,N3)  -X(1,N1))**2  +(X(2,N3)  -X(2,N1))**2 &
                                               +(X(3,N3)  -X(3,N1))**2)
      C1        = PROD1/(DA*DB)
      CHALF     = .5*(1.  +C1)
      IF (CHALF.LT.1.E-6) GO TO 200
      T1        = SQRT(1./CHALF  -1)
      ANGL1     = 2.*ATAN(T1)
      S1        = SIN(ANGL1)
      ANGL1     = RAD*ANGL1
      PROD2     = (X(1,N3)  -X(1,N2))*(X(1,N1)  -X(1,N2)) &
                 +(X(2,N3)  -X(2,N2))*(X(2,N1)  -X(2,N2)) &
                 +(X(3,N3)  -X(3,N2))*(X(3,N1)  -X(3,N2))
      DA        = SQRT((X(1,N3)  -X(1,N2))**2  +(X(2,N3)  -X(2,N2))**2 &
                                               +(X(3,N3)  -X(3,N2))**2)
      DB        = SQRT((X(1,N1)  -X(1,N2))**2  +(X(2,N1)  -X(2,N2))**2 &
                                               +(X(3,N1)  -X(3,N2))**2)
      C2        = PROD2/(DA*DB)
      CHALF     = .5*(1.  +C2)
      IF (CHALF.LT.1.E-6) GO TO 200
      T2        = SQRT(1./CHALF  -1)
      ANGL2     = 2.*ATAN(T2)
      S2        = SIN(ANGL2)
      ANGL2     = RAD*ANGL2
      PROD3     = (X(1,N1)  -X(1,N3))*(X(1,N2)  -X(1,N3)) &
                 +(X(2,N1)  -X(2,N3))*(X(2,N2)  -X(2,N3)) &
                 +(X(3,N1)  -X(3,N3))*(X(3,N2)  -X(3,N3))
      DA        = SQRT((X(1,N1)  -X(1,N3))**2  +(X(2,N1)  -X(2,N3))**2 &
                                               +(X(3,N1)  -X(3,N3))**2)
      DB        = SQRT((X(1,N2)  -X(1,N3))**2  +(X(2,N2)  -X(2,N3))**2 &
                                               +(X(3,N2)  -X(3,N3))**2)
      C3        = PROD3/(DA*DB)
      CHALF     = .5*(1.  +C3)
      IF (CHALF.LT.1.E-6) GO TO 200
      T3        = SQRT(1./CHALF  -1)
      ANGL3     = 2.*ATAN(T3)
      S3        = SIN(ANGL3)
      ANGL3     = RAD*ANGL3
      Q         = .5*(S1  +S2  +S3)/(S1*S2*S3)
      RETURN
  200 WRITE (6,600) J,N1,X(1,N1),X(2,N1),X(3,N1), &
                      N2,X(1,N2),X(2,N2),X(3,N2), &
                      N3,X(1,N3),X(2,N3),X(3,N3)
      STOP
  600 FORMAT(//5X,'FACE ',I6,' HAS AN ANGLE OF 180 DEGREES'// &
          5X,'VERTEX ',I6,'  X = ',F10.4,'  Y = ',F10.4,'  Z = ',F10.4/ &
          5X,'VERTEX ',I6,'  X = ',F10.4,'  Y = ',F10.4,'  Z = ',F10.4/ &
          5X,'VERTEX ',I6,'  X = ',F10.4,'  Y = ',F10.4,'  Z = ',F10.4// &
              5X,'PROGRAM STOPPED IN FANGLE')
      END SUBROUTINE FANGLE







!
!     ******************************************************************
!
      SUBROUTINE TETANG (X,NDC,NBH,IPROT,NCELL)
!
!     ******************************************************************
!     *                                                                *
!     *   CALCULATE MAXIMUM AND MINIMUM DIHEDRAL ANGLE AMONG ALL       *
!     *   TETRAHEDRA IN MESH                                           *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      INTEGER :: NCELL
      INTEGER :: NDC(4,*),NBH(4,*),IPROT(*)
      DOUBLE PRECISION :: X(3,*)
      
      INTEGER :: J,NPASS,N1,N2,N3,N4
      DOUBLE PRECISION :: ANGMAX,ANGMIN,ANG1,ANG2,ANG3,ANG4,ANG5,ANG6, &
                          DSMAX,DSMIN,DS1,DS2,DS3,DS4,DS5,DS6
      
!
!     ******************************************************************
!
      ANGMIN      = 180.
      ANGMAX      = 0.
      NPASS       = 0
      DO 20 J=1,NCELL
      IF (IPROT(J).EQ.1) GO TO 20
      IF (NDC(1,J).EQ.0) GO TO 20
      IF (NBH(1,J).EQ.0) GO TO 20
      N1          = NDC(1,J)
      N2          = NDC(2,J)
      N3          = NDC(3,J)
      N4          = NDC(4,J)
      IF (N4.EQ.-1) GO TO 20
      ANG1        = DIHED (N1,N2,N3,N4,X)
      IF (ANG1.GT.180.) ANG1 = 360.  -ANG1
      ANG2        = DIHED (N2,N3,N1,N4,X)
      IF (ANG2.GT.180.) ANG2 = 360.  -ANG2
      ANG3        = DIHED (N3,N1,N2,N4,X)
      IF (ANG3.GT.180.) ANG3 = 360.  -ANG3
      ANG4        = DIHED (N1,N4,N2,N3,X)
      IF (ANG4.GT.180.) ANG4 = 360.  -ANG4
      ANG5        = DIHED (N2,N4,N1,N3,X)
      IF (ANG5.GT.180.) ANG5 = 360.  -ANG5
      ANG6        = DIHED (N3,N4,N1,N2,X)
      IF (ANG6.GT.180.) ANG6 = 360.  -ANG6
      ANGMIN      = MIN(ANGMIN,ANG1,ANG2,ANG3,ANG4,ANG5,ANG6)
      ANGMAX      = MAX(ANGMAX,ANG1,ANG2,ANG3,ANG4,ANG5,ANG6)
      DS1         = (X(1,N1)  -X(1,N2))**2  +(X(2,N1) -X(2,N2))**2 &
                   +(X(3,N1)  -X(3,N2))**2
      DS2         = (X(1,N1)  -X(1,N3))**2  +(X(2,N1) -X(2,N3))**2 &
                   +(X(3,N1)  -X(3,N3))**2
      DS3         = (X(1,N1)  -X(1,N4))**2  +(X(2,N1) -X(2,N4))**2 &
                   +(X(3,N1)  -X(3,N4))**2
      DS4         = (X(1,N2)  -X(1,N3))**2  +(X(2,N2) -X(2,N3))**2 &
                   +(X(3,N2)  -X(3,N3))**2
      DS5         = (X(1,N2)  -X(1,N4))**2  +(X(2,N2) -X(2,N4))**2 &
                   +(X(3,N2)  -X(3,N4))**2
      DS6         = (X(1,N3)  -X(1,N4))**2  +(X(2,N3) -X(2,N4))**2 &
                   +(X(3,N3)  -X(3,N4))**2
      IF (NPASS.EQ.1) GO TO 10
      NPASS       = 1
      DSMIN       = MIN(DS1,DS2,DS3,DS4,DS5,DS6)
      DSMAX       = MAX(DS1,DS2,DS3,DS4,DS5,DS6)
      GO TO 20
   10 DSMIN       = MIN(DSMIN,DS1,DS2,DS3,DS4,DS5,DS6)
      DSMAX       = MAX(DSMAX,DS1,DS2,DS3,DS4,DS5,DS6)
   20 CONTINUE
      DSMIN       = SQRT(DSMIN)
      DSMAX       = SQRT(DSMAX)
      WRITE (6,600) ANGMIN,ANGMAX,DSMIN,DSMAX
      RETURN
  600 FORMAT(  5X,'*************************************************'/ &
               5X,'*                                               *'/ &
               5X,'*  MINIMUM DIHEDRAL ANGLE IS ',F6.2,' DEGREES     *'/ &
               5X,'*  MAXIMUM DIHEDRAL ANGLE IS ',F6.2,' DEGREES     *'/ &
               5X,'*                                               *'/ &
               5X,'*  MINIMUM EDGE DISTANCE IS ',F12.4,'        *'/ &
               5X,'*  MAXIMUM EDGE DISTANCE IS ',F12.4,'        *'/ &
               5X,'*                                               *'/ &
               5X,'*************************************************')
      END SUBROUTINE TETANG






!
!     ******************************************************************
!
      SUBROUTINE TETCOF (N1,N2,N3,N4,XPT,YPT,ZPT,V1,V2,V3,V4,VTET,X)
!
!     ******************************************************************
!     *                                                                *
!     *  DETERMINE VOLUME VTET OF TETRAHEDRON (N1,N2,N3,N4) AND THE    *
!     *  VOLUMES OF THE FOUR TETRAHEDRA ASSOCIATED WITH ITS FACES      *
!     *  AND THE POINT (XPT,YPT,ZPT).                                  *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE
      
      DOUBLE PRECISION :: VTET,V1,V2,V3,V4,XPT,YPT,ZPT
      DOUBLE PRECISION :: X(3,*)

      INTEGER :: N1,N2,N3,N4
      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,RNX,RNY,RNZ
!
!     ******************************************************************
!
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      V4      = ABS(RNX*(XPT  -X(1,N1))  +RNY*(YPT  -X(2,N1)) &
                                         +RNZ*(ZPT  -X(3,N1)))
      VTET    = ABS(RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                             +RNZ*(X(3,N4)  -X(3,N1)))
      RNX     = COFACT (X(2,N2),X(2,N3),X(2,N4),X(3,N2),X(3,N3),X(3,N4))
      RNY     = COFACT (X(3,N2),X(3,N3),X(3,N4),X(1,N2),X(1,N3),X(1,N4))
      RNZ     = COFACT (X(1,N2),X(1,N3),X(1,N4),X(2,N2),X(2,N3),X(2,N4))
      V1      = ABS(RNX*(XPT  -X(1,N2))  +RNY*(YPT  -X(2,N2)) &
                                         +RNZ*(ZPT  -X(3,N2)))
      RNX     = COFACT (X(2,N3),X(2,N4),X(2,N1),X(3,N3),X(3,N4),X(3,N1))
      RNY     = COFACT (X(3,N3),X(3,N4),X(3,N1),X(1,N3),X(1,N4),X(1,N1))
      RNZ     = COFACT (X(1,N3),X(1,N4),X(1,N1),X(2,N3),X(2,N4),X(2,N1))
      V2      = ABS(RNX*(XPT  -X(1,N3))  +RNY*(YPT  -X(2,N3)) &
                                         +RNZ*(ZPT  -X(3,N3)))
      RNX     = COFACT (X(2,N4),X(2,N1),X(2,N2),X(3,N4),X(3,N1),X(3,N2))
      RNY     = COFACT (X(3,N4),X(3,N1),X(3,N2),X(1,N4),X(1,N1),X(1,N2))
      RNZ     = COFACT (X(1,N4),X(1,N1),X(1,N2),X(2,N4),X(2,N1),X(2,N2))
      V3      = ABS(RNX*(XPT  -X(1,N4))  +RNY*(YPT  -X(2,N4)) &
                                         +RNZ*(ZPT  -X(3,N4)))
      RETURN
      END SUBROUTINE TETCOF







!
!     ******************************************************************
!
      SUBROUTINE VOLCOM (L,NPT,VDIFF,NCONT,X,NDC)
!
!     ******************************************************************
!     *                                                                *
!     *  COMPUTE VOLUMES OF TETRAHEDRA (N1,N2,N3,N4), (N1,N2,N3,NPT)   *
!     *  (N2,N3,N4,NPT), (N3,N4,N1,NPT) AND (N4,N1,N2,NPT).            *
!     *  POINT NPT LIES INSIDE TETRAHEDRON (N1,N2,N3,N4) IF NCONT = 0  *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *    COPYRIGHT  (C)  TIM BAKER   1994                            *
!     *                                                                *
!     ******************************************************************
!
      IMPLICIT NONE

      DOUBLE PRECISION :: X(3,*)
      INTEGER :: L,N,NCNT,NCONT,NEND,NPT,N1,N2,N3,N4
      INTEGER :: NDC(4,*)

      DOUBLE PRECISION :: A1,A2,A3,B1,B2,B3,PTOT,RNX,RNY,RNZ,SC,SP,VDIFF            
      DOUBLE PRECISION :: V(4),VCELL(4),PTEST(4)
!
!     ******************************************************************
!
      N1        = NDC(1,L)
      N2        = NDC(2,L)
      N3        = NDC(3,L)
      N4        = NDC(4,L)
!
!     COMPUTE VOLUMES OF COMPONENT TETRAHEDRA
!
      NEND      = N4
      NCNT      = 0
   10 NCNT      = NCNT  +1
      RNX     = COFACT (X(2,N1),X(2,N2),X(2,N3),X(3,N1),X(3,N2),X(3,N3))
      RNY     = COFACT (X(3,N1),X(3,N2),X(3,N3),X(1,N1),X(1,N2),X(1,N3))
      RNZ     = COFACT (X(1,N1),X(1,N2),X(1,N3),X(2,N1),X(2,N2),X(2,N3))
      V(NCNT) = RNX*(X(1,NPT)  -X(1,N1))  +RNY*(X(2,NPT)  -X(2,N1)) &
                                          +RNZ*(X(3,NPT)  -X(3,N1))
      VCELL(NCNT) = RNX*(X(1,N4)  -X(1,N1))  +RNY*(X(2,N4)  -X(2,N1)) &
                                             +RNZ*(X(3,N4)  -X(3,N1))
      SP      = SIGN(1.0D0,V(NCNT))
      SC      = SIGN(1.0D0,VCELL(NCNT))
      V(NCNT) = SP*V(NCNT)
      VCELL(NCNT) = SC*VCELL(NCNT)
      PTEST(NCNT) = SP*SC
      IF (N1.EQ.NEND) GO TO 20
      N       = N1
      N1      = N2
      N2      = N3
      N3      = N4
      N4      = N
      GO TO 10
!
!     COMPUTE DIFFERENCE BETWEEN SUM OF COMPONENT TETRAHEDRA VOLUMES AND
!     VOLUME OF TETRAHEDRON THAT MAY CONTAIN POINT NPT
!
   20 VDIFF   = V(1)+V(2)+V(3)+V(4) &
              -0.25D0*(VCELL(1)+VCELL(2)+VCELL(3)+VCELL(4))
      PTOT    = PTEST(1)+PTEST(2)+PTEST(3)+PTEST(4)
      NCONT   = 0
      IF (PTOT.GT.3.0D0) RETURN
      NCONT   = 1
      IF (PTOT.GT.-3.0D0) RETURN
      NCONT   = -1

      WRITE (6,900) N1,X(1,N1),X(2,N1),X(3,N1), &
                    N2,X(1,N2),X(2,N2),X(3,N2), &
                    N3,X(1,N3),X(2,N3),X(3,N3), &
                    N4,X(1,N4),X(2,N4),X(3,N4), &
                    NPT,X(1,NPT),X(2,NPT),X(3,NPT), &
                    V(1),V(2),V(3),V(4), &
                    VCELL(1),VCELL(2),VCELL(3),VCELL(4), &
                   PTEST(1),PTEST(2),PTEST(3),PTEST(4)
  900 FORMAT('N1 ',I5,' X = ',F10.4,' Y = ',F10.4,'Z = ',F10.4/ & 
             'N2 ',I5,' X = ',F10.4,' Y = ',F10.4,'Z = ',F10.4/ &
             'N3 ',I5,' X = ',F10.4,' Y = ',F10.4,'Z = ',F10.4/ &
             'N4 ',I5,' X = ',F10.4,' Y = ',F10.4,'Z = ',F10.4/ &
             'NPT ',I5,' X = ',F10.4,' Y = ',F10.4,'Z = ',F10.4/ &
             'V1 = ',E13.5,' V2 = ',E13.5,' V3 = ',E13.5,' V4 = ',E13.5/ &
             'VCELL = ',4(2X,E13.5)/ &
             'P1 = ',E13.5,' P2 = ',E13.5,' P3 = ',E13.5,' P4 = ',E13.5)

      RETURN
      END SUBROUTINE VOLCOM


! ******************************************************************************
!     Compute cofactor (extracted into own function)
! ******************************************************************************

      FUNCTION COFACT(A1,A2,A3,B1,B2,B3) 
      
        IMPLICIT NONE
      
        DOUBLE PRECISION :: COFACT      
      
        DOUBLE PRECISION, INTENT(IN) :: A1,A2,A3,B1,B2,B3
        
        COFACT = (A2 -A1)*(B3 -B1) -(A3 -A1)*(B2 -B1)
     
      END FUNCTION COFACT


      END MODULE RFLU_ModRepair3D
  
! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModRepair3D.F90,v $
!   Revision 1.4  2008/12/06 08:44:23  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:34  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2003/07/09 22:38:56  haselbac
!   Removed NREP (not used)
!
!   Revision 1.1  2002/10/05 19:16:17  haselbac
!   Initial revision - NOT TESTED
! 
! ******************************************************************************
  







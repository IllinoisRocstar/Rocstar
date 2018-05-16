      SUBROUTINE STRUCT (NNODE,NCELL,NDC,NHB,IMIN,NMIN,IMAX,NMAX,NINC)
C
C     ******************************************************************
C     *                                                                *
C     *  GIVEN A VOLUME MESH OF TETRAHEDRA CREATE THE DATA STRUCTURE   *
C     *  (NHB(K,L), K=1,4) WHICH CONTAINS THE ADDRESSES OF THE FOUR    *
C     *  TETRAHEDRA THAT ARE NEIGHBORS TO TETRAHEDRON L.               *  
C     *                                                                *
C     ******************************************************************
C     ******************************************************************
C     *                                                                *
C     *    COPYRIGHT  (C)  TIM BAKER   2000                            *
C     *                                                                *
C     ******************************************************************
C
C     ******************************************************************
C
      DIMENSION   NDC(4,*),NHB(4,*),
     .            IMIN(*),NMIN(*),IMAX(*),NMAX(*),NINC(*)
C
C     ******************************************************************
C
C     INITIALIZE ARRAYS IMIN AND IMAX
C
      DO 10 N=1,NNODE
      IMIN(N)   = 0
      IMAX(N)   = 0
   10 CONTINUE
C
C     INITIALIZE ARRAY NINC
C
      DO 15 L=1,NCELL
      NINC(L)   = 0
   15 CONTINUE
C
C     GENERATE LINKED LISTS OF TETRAHEDRA
C
      DO 30 L=1,NCELL
      L1        = NDC(1,L)
      L2        = NDC(2,L)
      L3        = NDC(3,L)
      L4        = NDC(4,L)
      N1        = MIN(L1,L2,L3,L4)
      N2        = MIN(MAX(L1,L2),MAX(L1,L3),MAX(L1,L4),
     .                MAX(L2,L3),MAX(L2,L4),MAX(L3,L4))
      N3        = MAX(MIN(L1,L2),MIN(L1,L3),MIN(L1,L4),
     .                MIN(L2,L3),MIN(L2,L4),MIN(L3,L4))
      N4        = MAX(L1,L2,L3,L4)
      NMIN(L)   = 0
      NMAX(L)   = 0
      I         = IMIN(N1)
      IF (I.EQ.0) THEN
         IMIN(N1)   = L
      ELSE
   20    NEXT      = NMIN(I)
         IF (NEXT.NE.0) THEN
            I         = NEXT
            GO TO 20
         ELSE
            NMIN(I)   = L
         END IF
      END IF
      I         = IMAX(N4)
      IF (I.EQ.0) THEN
         IMAX(N4)     = L
      ELSE
   25    NEXT      = NMAX(I)
         IF (NEXT.NE.0) THEN
            I         = NEXT
            GO TO 25
         ELSE
            NMAX(I)  = L
         END IF
      END IF
   30 CONTINUE
C
C     CREATE NHB ARRAY BY EXAMINING THE FACES OF EACH TETRAHEDRON
C
      DO 70 L=1,NCELL
      L1        = NDC(1,L)
      L2        = NDC(2,L)
      L3        = NDC(3,L)
      L4        = NDC(4,L)
      DO 70 K=1,4
      N1        = MIN(L1,L2,L3)
      N2        = MIN(MAX(L1,L2),MAX(L2,L3),MAX(L3,L1))
      N3        = MAX(L1,L2,L3)
      I         = IMIN(N1)
      IF (I.EQ.0) GO TO 45
   40 J1        = NDC(1,I)
      J2        = NDC(2,I)
      J3        = NDC(3,I)
      J4        = NDC(4,I)
      K1        = MIN(J1,J2,J3,J4)
      K2        = MIN(MAX(J1,J2),MAX(J1,J3),MAX(J1,J4),
     .                MAX(J2,J3),MAX(J2,J4),MAX(J3,J4))
      K3        = MAX(MIN(J1,J2),MIN(J1,J3),MIN(J1,J4),
     .                MIN(J2,J3),MIN(J2,J4),MIN(J3,J4))
      K4        = MAX(J1,J2,J3,J4)
      IF (I.NE.L) THEN
         IF (K2.EQ.N2.AND.K3.EQ.N3) GO TO 55
         IF (K2.EQ.N2.AND.K4.EQ.N3) GO TO 55
         IF (K3.EQ.N2.AND.K4.EQ.N3) GO TO 55
      END IF
      I         = NMIN(I)
      IF (I.NE.0) GO TO 40
   45 I         = IMAX(N3)
      IF (I.EQ.0) GO TO 55
   50 J1        = NDC(1,I)
      J2        = NDC(2,I)
      J3        = NDC(3,I)
      J4        = NDC(4,I)
      K1        = MIN(J1,J2,J3,J4)
      K2        = MIN(MAX(J1,J2),MAX(J1,J3),MAX(J1,J4),
     .                MAX(J2,J3),MAX(J2,J4),MAX(J3,J4))
      K3        = MAX(MIN(J1,J2),MIN(J1,J3),MIN(J1,J4),
     .                MIN(J2,J3),MIN(J2,J4),MIN(J3,J4))
      K4        = MAX(J1,J2,J3,J4)
      IF (I.NE.L) THEN
         IF (K1.EQ.N1.AND.K2.EQ.N2) GO TO 55
         IF (K1.EQ.N1.AND.K3.EQ.N2) GO TO 55
         IF (K2.EQ.N1.AND.K3.EQ.N2) GO TO 55
      END IF
      I         = NMAX(I)
      IF (I.NE.0) GO TO 50
   55 NINC(L)   = NINC(L)  +1
      NHB(NINC(L),L) = I
      LL        = L1
      L1        = L2
      L2        = L3
      L3        = L4
      L4        = LL
   70 CONTINUE
C
C     CHECK WHETHER ALL FACES HAVE BEEN FOUND
C
      DO 90 L=1,NCELL
      IF (NINC(L).NE.4) GO TO 300
   90 CONTINUE
      RETURN
  300 WRITE (6,600) L,NINC(L)
      STOP
  600 FORMAT(//5X,'FOR TETRAHEDRON ',I7,' NINC = ',I2/
     .         5X,'NINC SHOULD EQUAL 4 FOR ALL TETRAHEDRA.'/
     .         5X,'PROGRAM STOPPED IN ROUTINE STRUCT.')
      END

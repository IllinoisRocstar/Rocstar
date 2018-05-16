      SUBROUTINE STRUCT_tri (NNODE,NCELL,NDC,NHB,
     .                       IMIN,NMIN,IMAX,NMAX,NINC)
C
C     ******************************************************************
C     *                                                                *
C     *  GIVEN A MESH OF TRIANGLES CREATE THE DATA STRUCTURE           *
C     *  (NHB(K,L), K=1,3) WHICH CONTAINS THE ADDRESSES OF THE THREE   *
C     *  TRIANGLES THAT ARE NEIGHBORS TO TRIANGLE L.                   *  
C     *                                                                *
C     ******************************************************************
C     ******************************************************************
C     *                                                                *
C     *    COPYRIGHT  (C)  TIM BAKER   2001                            *
C     *                                                                *
C     ******************************************************************
C
C     ******************************************************************
C
      DIMENSION   NDC(3,*),NHB(3,*),
     .            IMIN(*),NMIN(*),IMAX(*),NMAX(*),NINC(*)
      DIMENSION   NSHAKE(3)
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
C     GENERATE LINKED LISTS OF TRIANGLES
C
      DO 30 L=1,NCELL
      L1        = NDC(1,L)
      L2        = NDC(2,L)
      L3        = NDC(3,L)
      N1        = MIN(L1,L2,L3)
      N2        = MIN(MAX(L1,L2),MAX(L1,L3),MAX(L2,L3))
      N3        = MAX(L1,L2,L3)
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
      I         = IMAX(N3)
      IF (I.EQ.0) THEN
         IMAX(N3)     = L
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
C     CREATE NHB ARRAY BY EXAMINING THE FACES OF EACH TRIANGLE
C
      DO 70 L=1,NCELL
      L1        = NDC(1,L)
      L2        = NDC(2,L)
      L3        = NDC(3,L)
      DO 70 K=1,3
      N1        = MIN(L1,L2)
      N2        = MAX(L1,L2)
      I         = IMIN(N1)
      IF (I.EQ.0) GO TO 45
   40 J1        = NDC(1,I)
      J2        = NDC(2,I)
      J3        = NDC(3,I)
      K1        = MIN(J1,J2,J3)
      K2        = MIN(MAX(J1,J2),MAX(J1,J3),MAX(J2,J3))
      K3        = MAX(J1,J2,J3)
      IF (I.NE.L) THEN
         IF (K2.EQ.N2.OR.K3.EQ.N2) GO TO 55
      END IF
      I         = NMIN(I)
      IF (I.NE.0) GO TO 40
   45 I         = IMAX(N2)
      IF (I.EQ.0) GO TO 55
   50 J1        = NDC(1,I)
      J2        = NDC(2,I)
      J3        = NDC(3,I)
      K1        = MIN(J1,J2,J3)
      K2        = MIN(MAX(J1,J2),MAX(J1,J3),MAX(J2,J3))
      K3        = MAX(J1,J2,J3)
      IF (I.NE.L) THEN
         IF (K1.EQ.N1.OR.K2.EQ.N1) GO TO 55
      END IF
      I         = NMAX(I)
      IF (I.NE.0) GO TO 50
   55 NINC(L)   = NINC(L)  +1
      NHB(NINC(L),L) = I
      LL        = L1
      L1        = L2
      L2        = L3
      L3        = LL
   70 CONTINUE
C
C     CHECK WHETHER ALL FACES HAVE BEEN FOUND
C
      DO 90 L=1,NCELL
      IF (NINC(L).NE.3) GO TO 300
   90 CONTINUE
C
C     REORDER NHB ARRAY SO THAT THE INDEX FOR A GIVEN NEIGHBOR
C     CORRESPONDS TO THE INDEX OF THE OPPOSITE VERTEX THAT APPEARS
C     IN THE NDC ARRAY.
C
      DO 115 L=1,NCELL
      NSHAKE(1) = 0
      NSHAKE(2) = 0
      NSHAKE(3) = 0
      DO 110 J=1,3
      JP1       = MOD(J,3)  +1
      JP2       = MOD(JP1,3)  +1
      N1        = NDC(JP1,L)
      N2        = NDC(JP2,L)
      NMINN     = MIN (N1,N2)
      NMAXX     = MAX (N1,N2)
      DO 100 K=1,3
      L2        = NHB(K,L)
      IF (L2.GT.0) THEN
         M1        = NDC(1,L2)
         M2        = NDC(2,L2)
         M3        = NDC(3,L2)
         MEND      = M3
   95    MMINN     = MIN (M1,M2)
         MMAXX     = MAX (M1,M2)
         IF (MMINN.EQ.NMINN.AND.MMAXX.EQ.NMAXX) GO TO 105
         IF (M1.EQ.MEND) GO TO 100
         M         = M1
         M1        = M2
         M2        = M3
         M3        = M
         GO TO 95
      ENDIF
  100 CONTINUE
      GO TO 110
C
C     NEIGHBORING TETRAHEDRON L2 HAS BEEN FOUND
C
  105 NSHAKE(J) = L2
  110 CONTINUE
      NHB(1,L)  = NSHAKE(1)
      NHB(2,L)  = NSHAKE(2)
      NHB(3,L)  = NSHAKE(3)
  115 CONTINUE
      RETURN
  300 WRITE (6,600) L,NINC(L)
      STOP
  600 FORMAT(//5X,'FOR TRIANGLE ',I7,' NINC = ',I2/
     .         5X,'NINC SHOULD EQUAL 3 FOR ALL TRIANGLES.'/
     .         5X,'PROGRAM STOPPED IN ROUTINE STRUCT-tri.')
      END

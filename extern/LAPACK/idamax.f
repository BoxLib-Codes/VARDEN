      INTEGER FUNCTION IDAMAX(N,DX,INCX)
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
c
c  Purpose
c  =======
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
c     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC DABS

      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
c
c        code for increment not equal to 1
c
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END


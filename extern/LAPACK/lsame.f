      LOGICAL FUNCTION LSAME(CA,CB)
C     TEST IF TWO CHARACTERS ARE ESSENTIALLY THE SAME.
C     THE CHARACTER CB IS ONE OF THE FORTRAN SET. 
C     (LOWER AND UPPER CASE LETTERS ARE EQUIVALENT.)
C     THIS IS A SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C
C     THIS SUBPROGRAM IS MACHINE-DEPENDENT.
C     VERSION FOR CDC SYSTEMS USING 6-12 BIT REPRESENTATIONS.
      CHARACTER CA(*)
      CHARACTER *1 CB
      INTEGER ICIRFX
      DATA ICIRFX/62/
C     SEE IF THE FIRST CHAR. IN STRING CA EQUALS STRING CB. 
      LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
      IF (LSAME) RETURN
C     THE CHARS. ARE NOT IDENTICAL.  NOW CHECK THEM FOR EQUIVALENCE.
C     LOOK FOR THE 'ESCAPE' CHARACTER, CIRCUMFLEX, FOLLOWED BY
C     THE LETTER.
      IVAL = ICHAR(CA(2))
      IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN 
          LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
      END IF
*
      RETURN
C     END 
C     LOGICAL FUNCTION LSAME(CA,CB)
C     TEST IF TWO CHARACTERS ARE ESSENTIALLY THE SAME.
C     THE CHARACTER CB IS ONE OF THE FORTRAN SET. 
C     (LOWER AND UPPER CASE LETTERS ARE EQUIVALENT.)
C     THIS IS A SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C
C     THIS SUBPROGRAM IS MACHINE-DEPENDENT.
C     VERSION FOR ANY ASCII MACHINE.
C     CHARACTER *1 CA
C     CHARACTER *1 CB
C     INTEGER IOFF
C     DATA IOFF/32/ 
C     SEE IF THE  CHAR. IN STRING CA EQUALS STRING CB.
C     LSAME = CA .EQ. CB
C     IF (LSAME) RETURN
C     THE CHARS. ARE NOT IDENTICAL.  NOW CHECK THEM FOR EQUIVALENCE.
C     ISHIFT = ICHAR(CA) - IOFF
C     IF (ISHIFT.GE.ICHAR('A') .AND. ISHIFT.LE.ICHAR('Z')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     RETURN
C     END 
C
C     LOGICAL FUNCTION LSAME(CA,CB)
C     TEST IF TWO CHARACTERS ARE ESSENTIALLY THE SAME.
C     THE CHARACTER CB IS ONE OF THE FORTRAN SET. 
C     (LOWER AND UPPER CASE LETTERS ARE EQUIVALENT.)
C     THIS IS A SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C
C     THIS SUBPROGRAM IS MACHINE-DEPENDENT.
C     VERSION FOR ANY EBCDIC MACHINE.
C     CHARACTER *1 CA
C     CHARACTER *1 CB
C     INTEGER IOFF
C     DATA IOFF/64/ 
C     SEE IF THE  CHAR. IN STRING CA EQUALS STRING CB.
C     LSAME = CA .EQ. CB
C     IF (LSAME) RETURN
C     THE CHARS. ARE NOT IDENTICAL.  NOW CHECK THEM FOR EQUIVALENCE.
C     ISHIFT = ICHAR(CA) + IOFF
C     IF (ISHIFT.GE.ICHAR('A') .AND. ISHIFT.LE.ICHAR('I')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     IF (ISHIFT.GE.ICHAR('J') .AND. ISHIFT.LE.ICHAR('R')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     IF (ISHIFT.GE.ICHAR('S') .AND. ISHIFT.LE.ICHAR('Z')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     RETURN
C     END 
      END 

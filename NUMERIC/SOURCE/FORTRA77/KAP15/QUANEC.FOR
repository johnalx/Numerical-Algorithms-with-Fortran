C[BA*)
C[KA{F 15}{Numerical Integration}{Numerical Integration}*)
C[FE{F 15.3}{Newton-Cotes Formulas}{Newton-Cotes Formulas}*)
C[LE*)
      SUBROUTINE  QUANEC  (A,B,N,NRP,FCT,VAL,IEG,ERREST,IERR)
C[IX{QUANEC}*)
C
C*****************************************************************
C                                                                *
C  This subroutine integrates a given function FCT(X) over the   *
C  interval [A,B] using the summed NEWTON-COTES formulas.        *
C[BE*)
C  N indicates how many subdivisions are used for the NC-formula.*
C  Thus the step size becomes H1=(B-A)/(N*VN) and for this we    *
C  obtain the approximate value QH1 for the integral. In case
C  N > 1 the process is repeated with step size
C  H2=(B-A)/((N-1)*VN). A different H2 may be chosen as well; the*
C  program can be adjusted accordingly. An error estimate is com-*
C  puted for QH1 by comparison with the second approximation QH2.*
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A    : left endpoint of the interval.                         *
C  B    : right endpoint of the interval.                        *
C  N    : number N of the sub-intervals in the summed quadrature *
C         formula.                                               *
C  NRP  : index for the specific NEWTON-COTES formula to be used:*
C          = 1: trapezoidal rule                                 *
C          = 2: SIMPSON rule                                     *
C          = 3: 3/8 - formula                                    *
C          = 4: 4/90 - rule                                      *
C          = 5: 5/288 - rule                                     *
C          = 6: 6/840 - rule                                     *
C          = 7: 7/17280 - rule                                   *
C  FCT  : name of the FUNCTION FCT(X) that we want to integrate. *
C         It has to be provided by the user in the following     *
C         format:                                                *
C              DOUBLE PRECISION FUNCTION FCT (X).                *
C         The function has to be defined as EXTERNAL in the      *
C         calling program.                                       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  VAL    : computed approximate value of the integral.          *
C  IEG    : error order of the method chosen.                    *
C  ERREST : error estimate for N > 1.                            *
C  IERR   : error parameter:                                     *
C           = 0: everything o.k.                                 *
C           = 1: false input parameter.                          *
C           = 2: A > B, VAL is given the appropriate sign        *
C           = 3: A = B: VAL = 0.0                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  date     : 08.19.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CFORML(7,8), FRML(7), Q(2), H(2), IRANGE(7)
      INTEGER   IEG
C
      DATA CFORML /
     1   1.0D0,  1.0D0,  1.0D0,  7.0D0, 19.0D0,  41.0D0,  751.0D0,
     2   1.0D0,  4.0D0,  3.0D0, 32.0D0, 75.0D0, 216.0D0, 3577.0D0,
     3   0.0D0,  1.0D0,  3.0D0, 12.0D0, 50.0D0,  27.0D0, 1323.0D0,
     4   0.0D0,  0.0D0,  1.0D0, 32.0D0, 50.0D0, 272.0D0, 2989.0D0,
     5   0.0D0,  0.0D0,  0.0D0,  7.0D0, 75.0D0,  27.0D0, 2989.0D0,
     6   0.0D0,  0.0D0,  0.0D0,  0.0D0, 19.0D0, 216.0D0, 1323.0D0,
     7   0.0D0,  0.0D0,  0.0D0,  0.0D0,  0.0D0,  41.0D0, 3577.0D0,
     8   0.0D0,  0.0D0,  0.0D0,  0.0D0,  0.0D0,   0.0D0,  751.0D0  /
C
      DATA FRML   /   2.0D0,     6.0D0, 8.0D0, 90.0D0, 288.0D0,
     1              840.0D0, 17280.0D0 /
      DATA IRANGE / 2, 4, 4, 6, 6, 8, 8 /
      DATA IFLAG  / 0 /
C
      IF (IFLAG .EQ. 0) THEN
         IFLAG = 1
         DO 10 K = 1, 8
            DO 50 I = 1, 7
               CFORML(I,K) = CFORML(I,K) / FRML(I) * DBLE(I)
 50         CONTINUE
 10      CONTINUE
      ENDIF
C
C*  test the input parameters
C
      IERR = 1
      IF (NRP .LE. 0 .OR. NRP .GE. 8. .OR. N .LE. 0)  R E T U R N
C
      IF (A .EQ. B) THEN
         IERR = 3
         VAL  = 0.0D0
         R E T U R N
      ENDIF
C
C*  check the interval endpoints
C
      IERR = 0
      IF (A .GT. B) THEN
         IERR = 2
         BL   = B
         BU   = A
         FC   = -1.0D0
      ELSE
         BL = A
         BU = B
         FC = 1.0D0
      ENDIF
C
      IEG = IRANGE(NRP)
      DO 20 I = 1, 2
         KEND = N + 1 - I
         H(I) = (BU - BL) / DBLE(KEND * NRP)
         Q(I) = 0.0D0
         DO 30 K = 1, KEND
            XN = BL + DBLE((K-1) * NRP) * H(I)
            DO 40 L = 1, NRP+1
               Q(I) = Q(I)+CFORML(NRP,L)*FCT(XN+DBLE(L-1)*H(I))
 40         CONTINUE
 30      CONTINUE
         Q(I) = Q(I) * H(I)
         IF (N .EQ. 1) THEN
            VAL = Q(1) * FC
            R E T U R N
         ENDIF
 20   CONTINUE
C
      ERREST   = DABS( (Q(1)-Q(2)) / ((H(2)/H(1))**IEG - 1.0D0) )
      VAL = Q(1) * FC
C
      R E T U R N
      END

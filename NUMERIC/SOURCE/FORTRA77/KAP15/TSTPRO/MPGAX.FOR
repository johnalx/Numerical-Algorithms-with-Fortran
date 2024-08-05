      PROGRAM TEST
C*************************************************************************
C                                                                        *
C     Test program for the  SUBROUTINE GAX                               *
C                                                                        *
C*************************************************************************
C                                                                        *
C     Given  :       F(X) = -100(X-0.5) E(-50(X-0.5)^2) -2 E(-2X)        *
C                                                                        *
C                                                          1.0           *
C     To find :      Approximation Q for the integral I  =  I   F(X) DX  *
C                                                          0.0           *
C                                                                        *
C     Exact solution : I = -0.864664716763387+ , correct in 15 decimals  *
C                                                                        *
C     We compute the approximation adaptively in double precision.       *
C     Error bound  EPS = 1.D-06                                          *
C                                                                        *
C     The test example was computed thus:                                *
C                                                                        *
C[ I WEIGHTS  I   ERROR   I # NODES   I   INTEGRAL   I IFEHL I METHOD I ]*
C[ -------------------------------------------------------------------- ]*
C[ ==================================================================== ]*
C[ I     0    I  .394E-07 I      16   I -.8646647168 I    0  I TRAPEZ.I ]*
C[ ==================================================================== ]*
C[ I     2    I  .929E-08 I      16   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     3    I  .422E-06 I       4   I -.8646647173 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     4    I  .180E-06 I       2   I -.8646647523 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     5    I  .505E-09 I       2   I -.8646647169 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     6    I  .959E-12 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     7    I  .233E-14 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     8    I  .111E-14 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     9    I  .600E-14 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    10    I  .211E-13 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    11    I  .101E-13 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    12    I  .282E-13 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    13    I  .540E-13 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    14    I  .477E-13 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    15    I  .520E-12 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    16    I  .461E-13 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    17    I  .386E-12 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    18    I  .145E-12 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    19    I  .917E-11 I       2   I -.8646647168 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I    20    I  .298E-10 I       2   I -.8646647167 I    0  I GAUSS  I ]*
C[ ==================================================================== ]*
C[ I     2    I  .129E-07 I      16   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I     4    I  .176E-06 I       4   I -.8646647170 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I     6    I  .126E-07 I       2   I -.8646647192 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I     8    I  .125E-10 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    10    I  .123E-13 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    12    I  .000E+00 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    14    I  .555E-15 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    16    I  .666E-15 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    18    I  .666E-15 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    20    I  .444E-15 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I    22    I  .666E-15 I       2   I -.8646647168 I    0  I CLEN   I ]*
C[ ==================================================================== ]*
C[ I     0    I  .258E-10 I       2   I -.8646647168 I    0  I ROMBE. I ]*
C[ ==================================================================== ]*
C[ I     2    I  .129E-07 I      16   I -.8646647168 I    0  I NEWTON I ]*
C[ ==================================================================== ]*
C[ I     3    I  .492E-08 I      16   I -.8646647168 I    0  I NEWTON I ]*
C[ ==================================================================== ]*
C[ I     4    I  .440E-06 I       4   I -.8646647162 I    0  I NEWTON I ]*
C[ ==================================================================== ]*
C[ I     5    I  .248E-06 I       4   I -.8646647165 I    0  I NEWTON I ]*
C[ ==================================================================== ]*
C[ I     6    I  .205E-06 I       2   I -.8646646764 I    0  I NEWTON I ]*
C[ ==================================================================== ]*
C[ I     7    I  .126E-06 I       2   I -.8646646919 I    0  I NEWTON I ]*
C[ ==================================================================== ]*
C                                                                        *
C     Computer            : IBM-PC/AT with coprozessor                   *
C     Compiler            : MICROSOFT FORTRAN Ver 5.0                    *
C                                                                        *
C*************************************************************************
C                                                                        *
C  Author      : Norbert Vogt                                              *
C  Date        : 2.8.1990                                                *
C  Source code : FORTRAN 77                                                *
C                                                                        *
C*************************************************************************
C
C     Declarations
C
      PARAMETER    (NM=100, N = 1)
C
C     Functions
C
      EXTERNAL     FCT, MACHPD
C
C     Variables
C
      DOUBLE PRECISION     INTVAL(1:100,1:2), XN(1:NM*4),AK(0:100, 2)
      DOUBLE PRECISION     QWERT, EXE, EPS
      INTEGER              LQM, NMAX, TP(1:NM), ST(1:NM*2),IANZ
      integer              nn
C
      CHARACTER*7  A1 ,A2, A3, A4, A5
C
      DATA A1 ,A2 ,A3, A4, A5
     +     /'TRAPEZ.', 'GAUSS', 'CLEN', 'ROMBE.', 'NEWTON'/
C
C     initialize
C
C                desired accuracy
C      EPS         = 1.D-06
C                maximal number of interval halvings
C      NMAX        = NM
C                Interval
C      INTVAL(1,1) = 0.0  = left endpoint
C      INTVAL(1,2) = 1.0  = right endpoint
C
C
C     Output
C
      WRITE (*,900)
      WRITE (*,902)
C
C     loop over all methods
C
      DO 100 LQM = 1, 5
C
         EPS         = 1.D-06
         NMAX        = NM
         INTVAL(1,1) = 0.0D+00
         INTVAL(1,2) = 1.0D+00
         IF (LQM .EQ. 2 .OR. LQM .EQ. 3 .OR. LQM .EQ. 5) THEN
C
C           use  Gauss quadrature
C
            IF ( LQM .EQ. 2) THEN
C
C              compute integral with degrees 2, ...,  20
C
               DO 111 IANZ = 2, 20
                  nn = n
                  CALL GAX (INTVAL, EPS, nn, FCT, NMAX, LQM, XN,
     +                      QWERT, EXE, TP, ST, IFEHL, AK, IANZ)
                  WRITE (*,901) IANZ,EXE,NMAX,QWERT,IFEHL,A2
                  WRITE (*,902)
                  EPS         = 1.D-06
                  NMAX        = NM
                  INTVAL(1,1) =0.0D+00
                  INTVAL(1,2) =1.0D+00
111            CONTINUE
            ELSE
C
C              use  Clenshaw-Curtis formulas
C
               IF (LQM .EQ. 3) THEN
C
C                 compute integral with weights 2, ...,  22
C
                  DO 112 IANZ = 2, 22, 2
                     nn = n
                     CALL GAX (INTVAL, EPS, nn, FCT, NMAX, LQM, XN,
     +                         QWERT, EXE, TP, ST, IFEHL, AK, IANZ)
                     WRITE (*,901)IANZ, EXE, NMAX, QWERT, IFEHL,A3
                     WRITE (*,902)
                     EPS         = 1.D-06
                     NMAX        = NM
                     INTVAL(1,1) = 0.0D+00
                     INTVAL(1,2) = 1.0D+00
 112              CONTINUE
               ELSE
C
C                 use summed  Newton-Cotes formulas
C
                  DO 113 IANZ = 2, 7
C
C                    compute integral with 2, ...,  7 sums
C
                     nn = n
                     CALL GAX (INTVAL, EPS, nn, FCT, NMAX, LQM, XN,
     +                         QWERT, EXE, TP, ST,IFEHL, AK, IANZ)
                     WRITE (*,901)IANZ, EXE, NMAX, QWERT, IFEHL, A5
                     WRITE (*,902)
                     EPS         = 1.D-06
                     NMAX        = NM
                     INTVAL(1,1) = 0.0D+00
                     INTVAL(1,2) = 1.0D+00
 113              CONTINUE
               ENDIF
            ENDIF
         ELSEIF (LQM .EQ. 1) THEN
C
C           use trapezoidal rule
C
            IANZ = 0
            nn   = n
            CALL GAX (INTVAL, EPS, nn, FCT, NMAX, LQM, XN, QWERT,
     +                EXE, TP, ST, IFEHL, AK, IANZ)
            WRITE (*,901)IANZ, EXE, NMAX, QWERT, IFEHL, A1
            WRITE (*,902)
         ELSEIF (LQM .EQ. 4) THEN
C
C           use  Romberg quadrature
C
            IANZ = 0
            nn   = n
            CALL GAX (INTVAL, EPS, nn, FCT, NMAX, LQM, XN, QWERT,
     +                EXE, TP, ST, IFEHL, AK, IANZ)
            WRITE (*,901)IANZ, EXE, NMAX, QWERT, IFEHL, A4
            WRITE (*,902)
         ENDIF
 100  CONTINUE
      STOP
C
C     Format statements
C
900   FORMAT (1X, 'C[ I WEIGHTS  I   ERROR   I # NODES   I',
     +            '   INTEGRAL   I IFEHL I METHOD I',T74, ']*', /,
     +        1X, 'C[ ',68('-'), T74, ']*')
 901  FORMAT (1X, 'C[ I', 4X, I2, 4X, 'I', 1X, E9.3, 1X, 'I', 3X, I5,
     +             3X, 'I', 1X, F12.10, 1X, 'I', 2X, I3, 2X, 'I', 1X,
     +             A7, 'I', T74, ']*')
 902  FORMAT (1X,'C[ ',68('='), T74, ']*')
      END
C
C
      DOUBLE PRECISION FUNCTION FCT(X)
C********************************************************************
C                                                                   *
C     Test function for SUBROUTINE GAX                              *
C                                                                   *
C     F(X) = -100(X - 0.5) E (-50(X-0.5)^2) -2 E (-2 * X)           *
C                                                                   *
C********************************************************************
      DOUBLE PRECISION X
      FCT = -1.D+02*(X-0.5D+00)*DEXP(-5.D+01*(X-0.5D+00)*(X-0.5D+00))
     +      -2.D+00* DEXP(-2.D+00*X)
      RETURN
      END

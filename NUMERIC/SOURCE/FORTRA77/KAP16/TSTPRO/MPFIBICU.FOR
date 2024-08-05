      PROGRAM TEST
C
C***************************************************************************
C                                                                          *
C   Test program for the subroutines BICSP2 and FIBICU      03/29/89       *
C                                                                          *
C     ( Elmar Pohl )                                                       *
C                                                                          *
C                                                                          *
C                                                                          *
C     Example for the use of  BICSP2 and FIBICU:                           *
C     We compute a bicubic spline surface with given normal vector and     *
C     integrate the spline function on the rectangle R = [XU,XO] x [YU,YO] *
C                                                                          *
C     (The output assumes fÅr M <= 10 )                                    *
C                                                                          *
C[                                                                        ]*
C[ COMPUTE A BICUBIC SPLINE USING  BICSP2                                 ]*
C[ ======================================                                 ]*
C[                                                                        ]*
C[ X-COORDINATES OF THE NODES                                             ]*
C[ --------------------------                                             ]*
C[                                                                        ]*
C[    I     X(I)                                                          ]*
C[ -----------------                                                      ]*
C[    0    0.000D+00                                                      ]*
C[    1    1.000D+00                                                      ]*
C[    2    2.000D+00                                                      ]*
C[    3    3.400D+00                                                      ]*
C[    4    5.000D+00                                                      ]*
C[                                                                        ]*
C[ Y-COORDINATES OF THE NODES                                             ]*
C[ --------------------------                                             ]*
C[                                                                        ]*
C[    I     Y(I)                                                          ]*
C[ -----------------                                                      ]*
C[    0   -1.500D+00                                                      ]*
C[    1    0.000D+00                                                      ]*
C[    2    1.000D+00                                                      ]*
C[    3    1.500D+00                                                      ]*
C[                                                                        ]*
C[ GIVEN FUNCTION VALUES A(I,K,0,0) = S( X(I),Y(K) )                      ]*
C[ -------------------------------------------------                      ]*
C[            0           1           2           3                       ]*
C[ ---------------------------------------------------------------------- ]*
C[    0   0.000D+00   0.000D+00   0.000D+00   0.000D+00                   ]*
C[    1   7.074D-02   1.000D+00   5.403D-01   7.074D-02                   ]*
C[    2   2.829D-01   4.000D+00   2.161D+00   2.829D-01                   ]*
C[    3   8.177D-01   1.156D+01   6.246D+00   8.177D-01                   ]*
C[    4   1.768D+00   2.500D+01   1.351D+01   1.768D+00                   ]*
C[                                                                        ]*
C[ EXAMPLE FOR APPLICATIONS OF FIBICU:                                    ]*
C[ ===================================                                    ]*
C[                                                                        ]*
C[ INTEGRATION OF S(X,Y) OVER THE RECTANGLE:                              ]*
C[     R = [ X(0),X(N) ] X [ Y(0),Y(M) ]                                  ]*
C[ RESULT:   I(S;R) =  8.26914D+01                                        ]*
C                                                                          *
C***************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
      PARAMETER(N=4,M=3,KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),Y(0:M),F(80)

C
C     FORMAT VARIABLES
C
      CHARACTER FO1555*50, FO1560*50

C
C set up mesh points
C
      DATA (X(I),I=0,N)  /0.0D0,1.0D0,2.0D0,3.4D0,5.0D0/
      DATA (Y(I),I=0,M)  /-1.5D0,0.0D0,1.0D0,1.5D0/

C
C compute function values at mesh points
C
      DO 10 I=0,N
         DO 10 K=0,M
            A(I,K,0,0) = S(X(I),Y(K))
10    CONTINUE
C
      WRITE(*,1520)
      WRITE(*,1530) (I,X(I),I=0,N)
      WRITE(*,1540) (I,Y(I),I=0,M)
      WRITE(*,1550)
      WRITE (FO1555, 1555) M + 1
      WRITE(*,FO1555) (I,I=0,M)
      WRITE(*,1556)
      DO 20 I=0,N
           WRITE (FO1560, 1560) M + 1
           WRITE(*,FO1560) I,(A(I,K,0,0),K=0,M)
20    CONTINUE
      CALL BICSP2(N,M,A,X,Y,F,IFEHL)
      IF(IFEHL.NE.0) THEN
        WRITE(*,1500)
        STOP
      ENDIF
C
C
C     The coefficients of the spline have been computed and printed.
C     Now we apply  FIBICU:
C     -----------------------------------------------------
C
C
      CALL FIBICU(N,M,A,X,Y,WERT)
      WRITE(*,1610) WERT
C
      STOP
1500  FORMAT(1X, 'C[', T76, ']*', /,
     +       1X, 'C[ *** ERROR IN BICSP2 ***', T76, ']*')
1520  FORMAT(1X, 'C[', T76, ']*', /,
     +       1X, 'C[ COMPUTE A BICUBIC SPLINE USING  BICSP2',
     +           T76, ']*', /,
     +       1X, 'C[ ', 38('='), T76, ']*')
1530  FORMAT(1X, 'C[', T76, ']*', /,
     +       1X, 'C[ X-COORDINATES OF THE NODES', T76, ']*', /,
     +       1X, 'C[ ', 26('-'), T76, ']*', /,
     +       1X, 'C[', T76, ']*', /,
     +       1X, 'C[    I', 5X, 'X(I)', T76, ']*', /,
     +       1X, 'C[ ', 17('-'), T76, ']*', /,
     +      (1X, 'C[ ', 1P, I4, 3X, D10.3, T76, ']*'))
1540  FORMAT(1X, 'C[', T76, ']*', /,
     +       1X, 'C[ Y-COORDINATES OF THE NODES', T76, ']*', /,
     +       1X, 'C[ ', 26('-'), T76, ']*', /,
     +       1X, 'C[', T76, ']*', /,
     +       1X, 'C[    I', 5X, 'Y(I)', T76, ']*', /,
     +       1X, 'C[ ', 17('-'), T76, ']*', /,
     +      (1X, 'C[ ', 1P, I4, 3X, D10.3, T76, ']*'))
1550  FORMAT(1X, 'C[', T76, ']*', /,
     +       1X, 'C[ GIVEN FUNCTION VALUES A(I,K,0,0) = S( X(I),',
     +           'Y(K) )', T76, ']*', /,
     +       1X, 'C[ ', 49('-'),  T76, ']*')
1555  FORMAT('(1X, ''C[ '',', I2, '(8X,I4), T76, '']*'')')
1556  FORMAT(1X, 'C[ ', 70('-'), T76, ']*')
1560  FORMAT('(1X, ''C[ '', I4, 1P, ', I2, '(2X,D10.3), T76, '']*'')')
1610  FORMAT(1X, 'C[', T76, ']*', /,
     +       1X, 'C[ EXAMPLE FOR APPLICATIONS OF FIBICU:',T76, ']*', /,
     +       1X, 'C[ ', 35('='), T76, ']*', /,
     +       1X, 'C[', T76, ']*', /,
     +       1X, 'C[ INTEGRATION OF S(X,Y) OVER THE RECTANGLE: ',
     +            T76, ']*', /,
     +       1X, 'C[     R = [ X(0),X(N) ] X [ Y(0),Y(M) ]', T76,']*',/,
     +       1X, 'C[ RESULT:   I(S;R) = ', 1PD12.5, T76, ']*')
      END
C
      DOUBLE PRECISION FUNCTION S (X,Y)
C
C*****************************************************************
C                                                                *
C  COMPUTE FUNCTIONAL VALUES
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION X,Y
      S = X*X*DCOS(Y)
      RETURN
      END

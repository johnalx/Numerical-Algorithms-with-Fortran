C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 12.2}
C[  {Two-Dimensional Interpolating Surface Splines}
C[  {Two-Dimensional Interpolating Surface Splines}*)
C[LE*)
      SUBROUTINE PROB2 (NX,X,Y,F,M,MARK,C,A,IWORK,WK)
C[IX{PROB2}*)
C
C*****************************************************************
C                                                                *
C PROB2 determines 2-dimensional surface splines for a given set *
C of triples (X(I),Y(I), F(X(I),Y(I)), I=1, ..., NX. The points  *
C (X(I),Y(I)) must be distinct. For each pair (X,Y) there must   *
C be precisely one value F=F(X,Y), i.e., F has to be a function. *
C The nodes do not have to lie on a rectangular grid and may be  *
C ordered in any way. It is advisable to transform the nodes     *
C (X(I),Y(I)) to the unit circle. A program for this is included.*
C[BE*)
C The derivative order should not be chosen too high, since the  *
C condition of the system of linear equations that must be       *
C solved worsens with increasing derivative order. Tests have    *
C shown that derivative orders between 3 and 5 are preferred.    *
C Only in a few of our tests have higher derivative orders       *
C resulted in marked improvements. The condition of the system   *
C matrix will also deteriorate with an increasing number of      *
C nodes or a decreasing distance between the nodes.              *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C NX   : number of nodes                                         *
C X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *
C        function values are known                               *
C F    : vector F(1:NX); the function values at the given nodes  *
C M    : derivative order for which the coefficients are to be   *
C        determined                                              *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C C    : vector C(1:(NX + M*(M+1)/2)); the coefficients of the   *
C        spline function                                         *
C MARK : indicates whether the system of equations is solvable   *
C        MARK = 1:  everything o.k.                              *
C        MARK = 0:  system matrix is numerically singular        *
C                                                                *
C                                                                *
C AUXILIARY VARIABLES:                                           *
C ====================                                           *
C A    : vector (1:(NX + M*(M+1)/2)*(3 + NX + M*(M+1)/2)/2)      *
C IWORK: INTEGER vector IWORK(1:(NX + M*(M+1)/2))                *
C WK   : vector WK(1:(NX+M*(M+1)/2)*((NX+M*(M+1)/2)+1)/2)        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C subroutines required: ALPHA2, GAMMA2, NEXT2, E2, CEPSPM,       *
C                       SESSPM                                   *
C                                                                *
C*****************************************************************
C                                                                *
C author   : Richard Reuter, 1983                                *
C editor   : Hartmut Turowski                                    *
C date     : 06.10.1988                                          *
C source   : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   declarations
C..
      DIMENSION  X(NX),Y(NX),F(NX),C(NX+M*(M+1)/2),
     F           A((NX+M*(M+1)/2)*(3+NX+M*(M+1)/2)/2),
     F           WK((NX+M*(M+1)/2)*((NX+M*(M+1)/2)+1)/2)
      INTEGER IWORK(NX+M*(M+1)/2)
C..
C..   order of the matrix
C..
      NM = NX + M*(M+1)/2
C..
C..   indicator for the polynomial components of the matrix
C..
      NXX = 1 + NX*(NX+1)/2
C..
C..   initializing the error parameter
C..
      MARK = 1
C..
C..   forming the matrix:
C..   polynomial components P: top right, in condensed form
C..
      CALL ALPHA2(NX,X,Y,M,A(NXX),IWORK(1),IWORK(1+(M+1)*M/2))
C..
C..   kernel component G of the matrix:
C..   top left, upper triangle in condensed form
C..
      CALL GAMMA2(NX,X,Y,M,A)
C..
C..   preparation of the right-hand side
C..
      DO 20 I = 1,NX
            C(I) = F(I)
   20 CONTINUE
      DO 30 I = NX+1, NM
            C(I) = 0.0D0
   30 CONTINUE
C..
C..   decomposition the system matrix of the linear equations
C..
      CALL CEPSPM(A,NM,C,IWORK,RCOND,A(NM*(NM+1)/2+1),WK)
C..
C..   if the matrix is numerically singular: stop
C..
      IF (1.0D0 .EQ. 1.0D0+RCOND) THEN
         MARK = 0
         RETURN
      ENDIF
C..
C..   solving the system of equations
C..
      CALL SESSPM(WK,NM,IWORK,C)
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE ALPHA2 (NX,X,Y,M,A,IDX,IDY)
C[IX{ALPHA2}*)
C
C*****************************************************************
C                                                                *
C Determine the polynomial components P of the syatem matrix.    *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C NX   : number of nodes                                         *
C X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *
C        function values are known                               *
C M    : derivative order for which the coefficients are to      *
C        be determined                                           *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C A    : vector                                                  *
C        A(1:(NX * (M*(M+1)/2) + (M*(M+1)/2 * (1+M*(M+1)/2)/2)));*
C        the polynomial components of the matrix in condensed    *
C        form                                                    *
C                                                                *
C                                                                *
C AUXILIARY VARIABLES:                                           *
C ====================                                           *
C IDX  : INTEGER vector IDX(1:((M+1)*(M+2)/2))                   *
C IDY  : INTEGER vector IDY(1:((M+1)*(M+2)/2))                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C subroutines required: NEXT2                                    *
C                                                                *
C*****************************************************************
C                                                                *
C author   : Richard Reuter, 1983                                *
C editor   : Hartmut Turowski                                    *
C date     : 06.10.1988                                          *
C source   : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   declarations
C..
      DIMENSION X(NX),Y(NX),
     F          A(NX*M*(M+1)/2+(M*(M+1)/2*(1+M*(M+1)/2)/2))
      INTEGER IDX((M+1)*(M+2)/2),IDY((M+1)*(M+2)/2)
C..
C..   the first monomial is 1.0
C..
      DO 10 I = 1, NX
         A(I) = 1.0D0
   10 CONTINUE
      A(NX+1) = 0.0D0
      IDX(1) = 0
      IDY(1) = 0
      DO 50 I = 2, M*(M+1)/2
C..
C..   determine the index of the monomial that is to be
C..   multiplied by X or Y
C..
         CALL NEXT2(I-1,IDX,IDY,IXY,K01)
         KL = IXY*(IXY-1)/2 + NX*(IXY-1)
         KLI= I  *(I  -1)/2 + NX*(I  -1)
         IF (K01 .EQ. 1) THEN
            DO 20 J = 1, NX
               A(KLI+J) = A(KL+J)*X(J)
   20       CONTINUE
         ELSE
            DO 30 J = 1, NX
               A(KLI+J) = A(KL+J)*Y(J)
   30       CONTINUE
         ENDIF
C..
C..    assign zero to the remaining matrix elements
C..
         DO 40 J = KLI+NX+1, KLI+NX+I
            A(J) = 0.0D0
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE NEXT2 (I,IDX,IDY,IXY,K)
C[IX{NEXT2}*)
C
C*****************************************************************
C                                                                *
C Auxiliary routine for efficiently determining all 2-dimensional*
C monomials up to degree M; see SUBROUTINE ALPHA2.               *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C I   : index of the 2-dimensional monomial determined previously*
C IDX : vector IDX(1:(I+1)); powers of X of the monomials with   *
C       index 1 to I                                             *
C IDY : vector IDY(1:(I+1)); powers of Y of the monomials with   *
C       index 1 to I                                             *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C IDX : vector IDX(1:(I+1)); IDX(I+1) is the power of X of the   *
C       monomial with index I+1                                  *
C IDY : vector IDY(1:(I+1)); IDY(I+1) is the power of Y of the   *
C       monomial with index I+1                                  *
C IXY : index of the monomial, that is to be multiplied by X or  *
C       Y in order to obtain the (I+1)-st monomial.              *
C K   : switching variable, multiplication by X or by Y          *
C       K=1 :  monom(I+1) = monom(IXY)*X                         *
C       K=0 :  monom(I+1) = monom(IXY)*Y                         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C subroutines required: none                                     *
C                                                                *
C*****************************************************************
C                                                                *
C author   : Richard Reuter, 1983                                *
C editor   : Hartmut Turowski                                    *
C date     : 06.10.1988                                          *
C source   : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C..
C..   declarations
C..
      INTEGER IDX(I+1),IDY(I+1)
C..
      N = IDX(I) + IDY(I)
      IF (IDX(I) .EQ. 0) THEN
         IDX(I+1) = N + 1
         IDY(I+1) = 0
         DO 10 J = 1, I
            IF (IDX(J) .EQ. N) GOTO 20
   10    CONTINUE
   20    CONTINUE
         IXY = J
         K   = 1
      ELSE
         IDX(I+1) = IDX(I) - 1
         IDY(I+1) = IDY(I) + 1
         DO 30 J = 1, I
            IF ((IDX(J) .EQ. IDX(I)-1) .AND.
     F          (IDY(J) .EQ. IDY(I))) GOTO 40
   30    CONTINUE
   40    CONTINUE
         IXY = J
         K   = 0
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GAMMA2 (NX,X,Y,M,A)
C[IX{GAMMA2}*)
C
C*****************************************************************
C                                                                *
C Initializing the kernel function components G of the matrix.   *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C NX   : number of nodes                                         *
C X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *
C        coefficients are being determined                       *
C M    : derivative order for which the coefficients are         *
C        determined                                              *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C A    : vector A(1:(NX*(NX+1)/2)); kernel function components   *
C        of the system matrix in condensed form                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C subroutines required: E2                                       *
C                                                                *
C*****************************************************************
C                                                                *
C author   : Richard Reuter, 1983                                *
C editor   : Hartmut Turowski                                    *
C date     : 06.10.1988                                          *
C source   : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   declarations
C..
      DIMENSION X(NX),Y(NX),A(NX*(NX+1)/2)
C..
C..   determine the kernel function components G
C..
      L = 0
      DO 20 I = 1, NX
         DO 10 K = 1, I - 1
            L = L + 1
            A(L) = E2(X(K)-X(I),Y(K)-Y(I),M)
   10    CONTINUE
         L = L + 1
C..
C..   set the main diagonal of G equal to zero;
C..   (we are dealing with interpolating splines)
C..
         A(L) = 0.0D0
   20 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION E2 (X,Y,M)
C[IX{E2}*)
C
C*****************************************************************
C                                                                *
C Evaluation of E at (X,Y). (See formula (12.19))                *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C X,Y : point where evaluation is to be performed                *
C M   : derivative order                                         *
C                                                                *
C                                                                *
C OUTPUT PARAMETER:                                              *
C =================                                              *
C E2  : functional value in DOUBLE PRECISION                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C subroutines required: none                                     *
C                                                                *
C*****************************************************************
C                                                                *
C author   : Richard Reuter                                      *
C editor   : Hartmut Turowski                                    *
C date     : 06.10.1988                                          *
C source   : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   determine the kernel function
C..
      R2 = X*X + Y*Y
      IF (R2 .EQ. 0.0D0) THEN
         E2 = 0.0D0
      ELSE
         E2 = DLOG(R2) * R2**(M-1)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE APPRX2 (X0,Y0,NX,M,X,Y,C,AP)
C[IX{APPRX2}*)
C
C*****************************************************************
C                                                                *
C Evaluation function for the interpolation.                     *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETER:                                               *
C ================                                               *
C X0,Y0: point where the function is to be evaluated             *
C NX   : number of nodes                                         *
C M    : derivative order for which the coefficients were        *
C        determined                                              *
C X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *
C        coefficients were determined                            *
C C    : vector C(1:(NX + M*(M+1)/2)); coefficient vector        *
C                                                                *
C                                                                *
C OUTPUT PARAMETER:                                              *
C =================                                              *
C AP   : approximate value for E at (X0,Y0)                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C subroutines required: none                                     *
C                                                                *
C*****************************************************************
C                                                                *
C author   : Richard Reuter, 1983                                *
C editor   : Hartmut Turowski                                    *
C date     : 06.10.1988                                          *
C source   : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   declarations
C..
      DIMENSION X(NX),Y(NX),C(NX+M*(M+1)/2)
C..
C..   for differing M the evaluations are performed differently:
C..   1. M = 1, 2, 3 ; specially coded, hence fast
C..   2. M > 3 ; each monomial is represented in the form
C..              (X**IX)*(Y**IY). The evaluation is slow and prone
C..              to rounding errors.
C..
C..   the first polynomial is always taken equal to 1
C..
      AP = C(NX+1)
      IF (M .EQ. 1) GOTO 20
      IF (M .EQ. 2) THEN
         AP = AP + C(NX+2)*X0 + C(NX+3)*Y0
      ELSE IF (M .EQ. 3) THEN
         AP = AP + (C(NX+2) + C(NX+4)*X0 + C(NX+5)*Y0)*X0
     F           + (C(NX+3) + C(NX+6)*Y0)*Y0
      ELSE
         IX = 0
         IY = 0
         DO 10 I = 2, M*(M+1)/2
            IF (IX .EQ. 0) THEN
               IX = IY + 1
               IY = 0
                  AP = AP + C(NX+I)*(X0**IX)
            ELSE
               IX = IX - 1
               IY = IY + 1
               IF (IX .EQ. 0) THEN
                     AP = AP + C(NX+I)*(Y0**IY)
               ELSE
                     AP = AP + C(NX+I)*(X0**IX)*(Y0**IY)
               ENDIF
            ENDIF
   10    CONTINUE
      ENDIF
   20 CONTINUE
C..
C..   component of the kernel function E
C..
C..   The function E2(X,Y,M) could be called at this point.
C..   However, this would slow down the evaluation considerably.
C..   Thus, a direct code is performed.
C..
      DO 30 I = 1, NX
            R2 = (X(I) - X0)**2 + (Y(I) - Y0)**2
            IF (R2 .EQ. 0.0D0) R2 = 1.0D0
            AP = AP + C(I)*DLOG(R2)*R2**(M-1)
   30 CONTINUE
      RETURN
      END

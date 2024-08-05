![          {Two--Dimensional Interpolating Surface Splines}*)          
      SUBROUTINE PROB2 (NX, X, Y, F, M, MARK, C, A, IWORK, WK) 
!                                                                       
!*****************************************************************      
!                                                                *      
! PROB2 determines 2-dimensional surface splines for a given set *      
! of triples (X(I),Y(I), F(X(I),Y(I)), I=1, ..., NX. The points  *      
! (X(I),Y(I)) must be distinct. For each pair (X,Y) there must   *      
! be precisely one value F=F(X,Y), i.e., F has to be a function. *      
! The nodes do not have to lie on a rectangular grid and may be  *      
! ordered in any way. It is advisable to transform the nodes     *      
! (X(I),Y(I)) to the unit circle. A program for this is included.*      
! The derivative order should not be chosen too high, since the  *      
! condition of the system of linear equations that must be       *      
! solved worsens with increasing derivative order. Tests have    *      
! shown that derivative orders between 3 and 5 are preferred.    *      
! Only in a few of our tests have higher derivative orders       *      
! resulted in marked improvements. The condition of the system   *      
! matrix will also deteriorate with an increasing number of      *      
! nodes or a decreasing distance between the nodes.              *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   : number of nodes                                         *      
! X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *      
!        function values are known                               *      
! F    : vector F(1:NX); the function values at the given nodes  *      
! M    : derivative order for which the coefficients are to be   *      
!        determined                                              *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! C    : vector C(1:(NX + M*(M+1)/2)); the coefficients of the   *      
!        spline function                                         *      
! MARK : indicates whether the system of equations is solvable   *      
!        MARK = 1:  everything o.k.                              *      
!        MARK = 0:  system matrix is numerically singular        *      
!                                                                *      
!                                                                *      
! AUXILIARY VARIABLES:                                           *      
! ====================                                           *      
! A    : vector (1:(NX + M*(M+1)/2)*(3 + NX + M*(M+1)/2)/2)      *      
! IWORK: INTEGER vector IWORK(1:(NX + M*(M+1)/2))                *      
! WK   : vector WK(1:(NX+M*(M+1)/2)*((NX+M*(M+1)/2)+1)/2)        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: ALPHA2, GAMMA2, NEXT2, E2, CEPSPM,       *      
!                       SESSPM                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX), F (NX), C (NX + M * (M + 1) / 2),       &
      A ( (NX + M * (M + 1) / 2) * (3 + NX + M * (M + 1) / 2) / 2),     &
      WK ( (NX + M * (M + 1) / 2) * ( (NX + M * (M + 1) / 2) + 1)       &
      / 2)                                                              
      INTEGER IWORK (NX + M * (M + 1) / 2) 
!..                                                                     
!..   order of the matrix                                               
!..                                                                     
      NM = NX + M * (M + 1) / 2 
!..                                                                     
!..   indicator for the polynomial components of the matrix             
!..                                                                     
      NXX = 1 + NX * (NX + 1) / 2 
!..                                                                     
!..   initializing the error parameter                                  
!..                                                                     
      MARK = 1 
!..                                                                     
!..   forming the matrix:                                               
!..   polynomial components P: top right, in condensed form             
!..                                                                     
      CALL ALPHA2 (NX, X, Y, M, A (NXX), IWORK (1), IWORK (1 + (M + 1)  &
      * M / 2) )                                                        
!..                                                                     
!..   kernel component G of the matrix:                                 
!..   top left, upper triangle in condensed form                        
!..                                                                     
      CALL GAMMA2 (NX, X, Y, M, A) 
!..                                                                     
!..   preparation of the right-hand side                                
!..                                                                     
      DO 20 I = 1, NX 
         C (I) = F (I) 
   20 END DO 
      DO 30 I = NX + 1, NM 
         C (I) = 0.0D0 
   30 END DO 
!..                                                                     
!..   decomposition the system matrix of the linear equations           
!..                                                                     
      CALL CEPSPM (A, NM, C, IWORK, RCOND, A (NM * (NM + 1) / 2 + 1),   &
      WK)                                                               
!..                                                                     
!..   if the matrix is numerically singular: stop                       
!..                                                                     
      IF (1.0D0.EQ.1.0D0 + RCOND) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!..                                                                     
!..   solving the system of equations                                   
!..                                                                     
      CALL SESSPM (WK, NM, IWORK, C) 
      RETURN 
      END SUBROUTINE PROB2                          
!                                                                       
!                                                                       
      SUBROUTINE ALPHA2 (NX, X, Y, M, A, IDX, IDY) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Determine the polynomial components P of the syatem matrix.    *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   : number of nodes                                         *      
! X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *      
!        function values are known                               *      
! M    : derivative order for which the coefficients are to      *      
!        be determined                                           *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! A    : vector                                                  *      
!        A(1:(NX * (M*(M+1)/2) + (M*(M+1)/2 * (1+M*(M+1)/2)/2)));*      
!        the polynomial components of the matrix in condensed    *      
!        form                                                    *      
!                                                                *      
!                                                                *      
! AUXILIARY VARIABLES:                                           *      
! ====================                                           *      
! IDX  : INTEGER vector IDX(1:((M+1)*(M+2)/2))                   *      
! IDY  : INTEGER vector IDY(1:((M+1)*(M+2)/2))                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: NEXT2                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX), A (NX * M * (M + 1) / 2 + (M * (M + 1)  &
      / 2 * (1 + M * (M + 1) / 2) / 2) )                                
      INTEGER IDX ( (M + 1) * (M + 2) / 2), IDY ( (M + 1) * (M + 2)     &
      / 2)                                                              
!..                                                                     
!..   the first monomial is 1.0                                         
!..                                                                     
      DO 10 I = 1, NX 
         A (I) = 1.0D0 
   10 END DO 
      A (NX + 1) = 0.0D0 
      IDX (1) = 0 
      IDY (1) = 0 
      DO 50 I = 2, M * (M + 1) / 2 
!..                                                                     
!..   determine the index of the monomial that is to be                 
!..   multiplied by X or Y                                              
!..                                                                     
         CALL NEXT2 (I - 1, IDX, IDY, IXY, K01) 
         KL = IXY * (IXY - 1) / 2 + NX * (IXY - 1) 
         KLI = I * (I - 1) / 2 + NX * (I - 1) 
         IF (K01.EQ.1) THEN 
            DO 20 J = 1, NX 
               A (KLI + J) = A (KL + J) * X (J) 
   20       END DO 
         ELSE 
            DO 30 J = 1, NX 
               A (KLI + J) = A (KL + J) * Y (J) 
   30       END DO 
         ENDIF 
!..                                                                     
!..    assign zero to the remaining matrix elements                     
!..                                                                     
         DO 40 J = KLI + NX + 1, KLI + NX + I 
            A (J) = 0.0D0 
   40    END DO 
   50 END DO 
      RETURN 
      END SUBROUTINE ALPHA2                         
!                                                                       
!                                                                       
      SUBROUTINE NEXT2 (I, IDX, IDY, IXY, K) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Auxiliary routine for efficiently determining all 2-dimensional*      
! monomials up to degree M; see SUBROUTINE ALPHA2.               *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! I   : index of the 2-dimensional monomial determined previously*      
! IDX : vector IDX(1:(I+1)); powers of X of the monomials with   *      
!       index 1 to I                                             *      
! IDY : vector IDY(1:(I+1)); powers of Y of the monomials with   *      
!       index 1 to I                                             *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! IDX : vector IDX(1:(I+1)); IDX(I+1) is the power of X of the   *      
!       monomial with index I+1                                  *      
! IDY : vector IDY(1:(I+1)); IDY(I+1) is the power of Y of the   *      
!       monomial with index I+1                                  *      
! IXY : index of the monomial, that is to be multiplied by X or  *      
!       Y in order to obtain the (I+1)-st monomial.              *      
! K   : switching variable, multiplication by X or by Y          *      
!       K=1 :  monom(I+1) = monom(IXY)*X                         *      
!       K=0 :  monom(I+1) = monom(IXY)*Y                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: none                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
!..                                                                     
!..   declarations                                                      
!..                                                                     
      INTEGER IDX (I + 1), IDY (I + 1) 
!..                                                                     
      N = IDX (I) + IDY (I) 
      IF (IDX (I) .EQ.0) THEN 
         IDX (I + 1) = N + 1 
         IDY (I + 1) = 0 
         DO 10 J = 1, I 
            IF (IDX (J) .EQ.N) GOTO 20 
   10    END DO 
   20    CONTINUE 
         IXY = J 
         K = 1 
      ELSE 
         IDX (I + 1) = IDX (I) - 1 
         IDY (I + 1) = IDY (I) + 1 
         DO 30 J = 1, I 
            IF ( (IDX (J) .EQ.IDX (I) - 1) .AND. (IDY (J) .EQ.IDY (I) ) &
            ) GOTO 40                                                   
   30    END DO 
   40    CONTINUE 
         IXY = J 
         K = 0 
      ENDIF 
      RETURN 
      END SUBROUTINE NEXT2                          
!                                                                       
!                                                                       
      SUBROUTINE GAMMA2 (NX, X, Y, M, A) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Initializing the kernel function components G of the matrix.     *    
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   : number of nodes                                         *      
! X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *      
!        coefficients are being determined                       *      
! M    : derivative order for which the coefficients are         *      
!        determined                                              *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! A    : vector A(1:(NX*(NX+1)/2)); kernel function components   *      
!        of the system matrix in condensed form                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: E2                                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX), A (NX * (NX + 1) / 2) 
!..                                                                     
!..   determine the kernel function components G                        
!..                                                                     
      L = 0 
      DO 20 I = 1, NX 
         DO 10 K = 1, I - 1 
            L = L + 1 
            A (L) = E2 (X (K) - X (I), Y (K) - Y (I), M) 
   10    END DO 
         L = L + 1 
!..                                                                     
!..   set the main diagonal of G equal to zero;                         
!..   (we are dealing with interpolating splines)                       
!..                                                                     
         A (L) = 0.0D0 
   20 END DO 
      RETURN 
      END SUBROUTINE GAMMA2                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION E2 (X, Y, M) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Evaluation of E at (X,Y). (See formula (12.19))                *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X,Y : point where evaluation is to be performed                *      
! M   : derivative order                                         *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! E2  : functional value in DOUBLE PRECISION                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: none                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter                                      *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   determine the kernel function                                     
!..                                                                     
      R2 = X * X + Y * Y 
      IF (R2.EQ.0.0D0) THEN 
         E2 = 0.0D0 
      ELSE 
         E2 = DLOG (R2) * R2** (M - 1) 
      ENDIF 
      RETURN 
      END FUNCTION E2                               
!                                                                       
!                                                                       
      SUBROUTINE APPRX2 (X0, Y0, NX, M, X, Y, C, AP) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Evaluation function for the interpolation.                     *      
!                                                                *      
!                                                                *      
! INPUT PARAMETER:                                               *      
! ================                                               *      
! X0,Y0: point where the function is to be evaluated             *      
! NX   : number of nodes                                         *      
! M    : derivative order for which the coefficients were        *      
!        determined                                              *      
! X,Y  : vectors X(1:NX), Y(1:NX); the nodes for which the       *      
!        coefficients were determined                            *      
! C    : vector C(1:(NX + M*(M+1)/2)); coefficient vector        *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! AP   : approximate value for E at (X0,Y0)                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: none                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX), C (NX + M * (M + 1) / 2) 
!..                                                                     
!..   for differing M the evaluations are performed differently:        
!..   1. M = 1, 2, 3 ; specially coded, hence fast                      
!..   2. M > 3 ; each monomial is represented in the form               
!..              (X**IX)*(Y**IY). The evaluation is slow and prone      
!..              to rounding errors.                                    
!..                                                                     
!..   the first polynomial is always taken equal to 1                   
!..                                                                     
      AP = C (NX + 1) 
      IF (M.EQ.1) GOTO 20 
      IF (M.EQ.2) THEN 
         AP = AP + C (NX + 2) * X0 + C (NX + 3) * Y0 
      ELSEIF (M.EQ.3) THEN 
         AP = AP + (C (NX + 2) + C (NX + 4) * X0 + C (NX + 5) * Y0)     &
         * X0 + (C (NX + 3) + C (NX + 6) * Y0) * Y0                     
      ELSE 
         IX = 0 
         IY = 0 
         DO 10 I = 2, M * (M + 1) / 2 
            IF (IX.EQ.0) THEN 
               IX = IY + 1 
               IY = 0 
               AP = AP + C (NX + I) * (X0**IX) 
            ELSE 
               IX = IX - 1 
               IY = IY + 1 
               IF (IX.EQ.0) THEN 
                  AP = AP + C (NX + I) * (Y0**IY) 
               ELSE 
                  AP = AP + C (NX + I) * (X0**IX) * (Y0**IY) 
               ENDIF 
            ENDIF 
   10    END DO 
      ENDIF 
   20 CONTINUE 
!..                                                                     
!..   component of the kernel function E                                
!..                                                                     
!..   The function E2(X,Y,M) could be called at this point.             
!..   However, this would slow down the evaluation considerably.        
!..   Thus, a direct code is performed.                                 
!..                                                                     
      DO 30 I = 1, NX 
         R2 = (X (I) - X0) **2 + (Y (I) - Y0) **2 
         IF (R2.EQ.0.0D0) R2 = 1.0D0 
         AP = AP + C (I) * DLOG (R2) * R2** (M - 1) 
   30 END DO 
      RETURN 
      END SUBROUTINE APPRX2                         

      SUBROUTINE PROB3 (NX, X, Y, Z, F, M, MARK, C, A, IWORK, WK) 
!                                                                       
!*****************************************************************      
!                                                                *      
! PROB3 computes three-dimensional surface splines for arbitrary *      
! given points (X(I),Y(I),Z(I),F(X(I),Y(I)), I=1, ..., NX.       *      
! The nodes (X(I),Y(I),Z(I)) must be distinct and F must be a    *      
! function, i.e., for each (X,Y,Z) in the node set there must    *      
! correspond a unique  F=F(X,Y,Z). The nodes need not be ordered.*      
! The desired smoothness of the spline, i. e., its derivative    *      
! order should be stipulated as rather low since the condition   *      
! number of the system of equations that has to be solved worsens*      
! with increasing derivative order.                              *      
! Tests indicate that derivative orders between  3 and 5  can be *      
! recommended. Higher orders showed improvement only in rare     *      
! cases. For an increasing number of nodes, i. e., a decreasing  *      
! distance between the nodes the condition number of the linear  *      
! system also tends to worsen.                                   *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX     :  Number of nodes                                      *      
! X,Y,Z  :  NX-vectors ..(1:NX); the coordinates of the nodes    *      
! F      :  NX-vector F(1:NX); the functional values at the nodes*      
! M      :  derivative order used to determine the coefficients  *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
!                                                                *      
! C      :  vector C(1:(NX + M*(M+1)*(M+2)/6)); the coefficients *      
!           of the spline                                        *      
!                                                                *      
! MARK   :  indicates whether the linear system could be solved: *      
!           MARK = 1:  all is ok                                 *      
!           MARK = 0:  system matrix is numerically singular     *      
!                                                                *      
! AUXILIARY PARAMETERS:                                          *      
! =====================                                          *      
! A      : vector A(1:((NX + M*(M+1)*(M+2)/6) *                  *      
!                         * (NX + M*(M+1)*(M+2)/6 + 3))/2)       *      
! IWORK  : integer vector IWORK(1:(NX + M*(M+1)*(M+2)/6))        *      
! WK     : vector WK(1:((NX + M*(M+1)*(M+2)/6)*                  *      
!                         * (NX + M*(M+1)*(M+2)/6 + 1)/2))       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines: ALPHA3, GAMMA3, NEXT3, E3, CEPSPM,       *      
!                       ZSPMMK, PCOSOL, PCOLTG, SESSPM, SCAPRO,  *      
!                       VECMWC, ABSSUM, INDMAX, VECADD, VECXCH,  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (NX), Y (NX), Z (NX), F (NX) 
      DIMENSION C (NX + M * (M + 1) * (M + 2) / 6) 
      DIMENSION A ( (NX + M * (M + 1) * (M + 2) / 6) * (NX + M *        &
      (M + 1) * (M + 2) / 6 + 3) / 2)                                   
      DIMENSION WK ( (NX + M * (M + 1) * (M + 2) / 6) * ( (NX + M *     &
      (M + 1) * (M + 2) / 6) + 1) / 2)                                  
      DIMENSION IWORK (NX + M * (M + 1) * (M + 2) / 6) 
!..                                                                     
!..   Order of the matrix                                               
!..                                                                     
      M3 = M * (M + 1) * (M + 2) / 6 
      NM = NX + M3 
!..                                                                     
!..   Pointer for the polynomial part of the system matrix              
!..                                                                     
      NXX = 1 + (NX * (NX + 1) ) / 2 
!..                                                                     
!..   Initialize error parameter                                        
!..                                                                     
      MARK = 1 
!..                                                                     
!..   Form system matrix:                                               
!..   Polynomial part P appears in condensed form in upper right corner 
!..                                                                     
      CALL ALPHA3 (NX, X, Y, Z, M - 1, A (NXX), IWORK (1), IWORK (1 +   &
      M3), IWORK (1 + 2 * M3) )                                         
!..                                                                     
!..   G part of the matrix:                                             
!..   condensed in upper left corner                                    
!..                                                                     
      CALL GAMMA3 (NX, X, Y, Z, M, A) 
!..                                                                     
!..   Set up right hand side                                            
!..                                                                     
      DO 10 I = 1, NX 
         C (I) = F (I) 
   10 END DO 
      DO 20 I = NX + 1, NM 
         C (I) = 0.0D0 
   20 END DO 
!..                                                                     
!..   factor the system matrix                                          
!..                                                                     
      CALL CEPSPM (A, NM, C, IWORK, RCOND, A ( (NM * (NM + 1) ) / 2 + 1)&
      , WK)                                                             
!..                                                                     
!..   Stop if the system matrix is numerically singular                 
!..                                                                     
      IF (1.0D0.EQ.1.0D0 + RCOND) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!..                                                                     
!..   Solve the linear system                                           
!..                                                                     
      CALL SESSPM (WK, NM, IWORK, C) 
      RETURN 
      END SUBROUTINE PROB3                          
!                                                                       
!                                                                       
      SUBROUTINE ALPHA3 (NX, X, Y, Z, M, A, IDX, IDY, IDZ) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Computes the polynomial part P of the system matrix.           *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   : number of nodes                                         *      
! X,Y,Z: NX-vectors ..(1:NX); the coordinates of the nodes       *      
! M    : derivative order                                        *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! A    : polynomial part of the system matrix in condensed form: *      
!        A(1:(((M+1)*(M+2)*(M+3)/6) *                            *      
!                          * (2*NX+((M+1)*(M+2)*(M+3)/6)+1)/2)   *      
!                                                                *      
! AUXILIARY PARAMETERS:                                          *      
! =====================                                          *      
! IDX  : ]                                                       *      
! IDY  : ] vectors ..(1:((M+1)*(M+2)*(M+3)/6))                   *      
! IDZ  : ]                                                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines: NEXT3                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (NX), Y (NX), Z (NX) 
      DIMENSION A ( (M + 1) * (M + 2) * (M + 3) / 6 * (2 * NX + (M + 1) &
      * (M + 2) * (M + 3) / 6 + 1) / 2)                                 
      DIMENSION IDX ( ( (M + 1) * (M + 2) * (M + 3) / 6) ) 
      DIMENSION IDY ( ( (M + 1) * (M + 2) * (M + 3) / 6) ) 
      DIMENSION IDZ ( ( (M + 1) * (M + 2) * (M + 3) / 6) ) 
!..                                                                     
!..   the first monomial is  1.0                                        
!..                                                                     
      DO 10 I = 1, NX 
         A (I) = 1.0D0 
   10 END DO 
      A (NX + 1) = 0.0D0 
      IDX (1) = 0 
      IDY (1) = 0 
      IDZ (1) = 0 
      L = 1 
      DO 40 I = 1, M 
         DO 30 IX = I, 0, - 1 
            DO 20 IY = I - IX, 0, - 1 
               IZ = I - IY - IX 
               L = L + 1 
               IDX (L) = IX 
               IDY (L) = IY 
               IDZ (L) = IZ 
   20       END DO 
   30    END DO 
   40 END DO 
      DO 90 I = 2, L 
!..                                                                     
!..   determine the index of the monomial that needs                    
!..   to be multiplied by  X,Y or Z                                     
!..                                                                     
         CALL NEXT3 (I, IDX, IDY, IDZ, ID, K) 
         KL = ID * (ID-1) / 2 + NX * (ID-1) 
         KLI = I * (I - 1) / 2 + NX * (I - 1) 
         IF (K.EQ.1) THEN 
            DO 50 J = 1, NX 
               A (KLI + J) = A (KL + J) * X (J) 
   50       END DO 
         ELSEIF (K.EQ.2) THEN 
            DO 60 J = 1, NX 
               A (KLI + J) = A (KL + J) * Y (J) 
   60       END DO 
         ELSE 
            DO 70 J = 1, NX 
               A (KLI + J) = A (KL + J) * Z (J) 
   70       END DO 
         ENDIF 
!..                                                                     
!..   zero the rest of the system matrix                                
!..                                                                     
         DO 80 J = KLI + NX + 1, KLI + NX + I 
            A (J) = 0.0D0 
   80    END DO 
   90 END DO 
      RETURN 
      END SUBROUTINE ALPHA3                         
!                                                                       
!                                                                       
      SUBROUTINE NEXT3 (I, IDX, IDY, IDZ, ID, K) 
!                                                                       
!*****************************************************************      
!                                                                *      
! SUBROUTINE that efficiently determines all three-dimensional   *      
! monomials up to degree M; refer to SUBROUTINE ALPHA3           *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! I   :  Index of the two-dimensional monomial that was last     *      
!        computed                                                *      
! IDX :  ]  I-vectors ID..(1:I); the powers of X, Y or Z in the  *      
! IDY :  ]  monomials with index 1 to I                          *      
! IDZ :  ]                                                       *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! ID  :  Index of the monomial, that must be multiplied by X,Y   *      
!        or Z in order to obtain the I-th monomial               *      
! K   :  Switch that toggles multiplication by X,Y or Z          *      
!        K=1 : Monom(I) = Monom(ID)*X                            *      
!        K=2 : Monom(I) = Monom(ID)*Y                            *      
!        K=3 : Monom(I) = Monom(ID)*Z                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines:  none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION IDX (I), IDY (I), IDZ (I) 
!..                                                                     
      N = IDX (I) + IDY (I) + IDZ (I) 
      IF (IDX (I) .NE.0) THEN 
         K = 1 
         ID = I - (N * (N + 1) ) / 2 
      ELSEIF (IDY (I) .NE.0) THEN 
         K = 2 
         ID = I + 1 - ( (N + 1) * (N + 2) ) / 2 
      ELSE 
         K = 3 
         ID = I - ( (N + 1) * (N + 2) ) / 2 
      ENDIF 
      RETURN 
      END SUBROUTINE NEXT3                          
!                                                                       
!                                                                       
!                                                                       
      SUBROUTINE GAMMA3 (NX, X, Y, Z, M, A) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Initialize the G part of the system matrix.                    *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   :  number of nodes                                        *      
! X,Y,Z:  NX-vectors ..(1:NX); the coordinates of the nodes      *      
! M    :  derivative order                                       *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! A    :  vector A(1:(NX*(NX+1)/2)); the G part of the system    *      
!         matrix in condensed form                               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines: E2                                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (NX), Y (NX), Z (NX), A (NX * (NX + 1) / 2) 
!..                                                                     
!..   Calculate the G part                                              
!..                                                                     
      L = 0 
      DO 20 I = 1, NX 
         DO 10 K = 1, I - 1 
            L = L + 1 
!..                                                                     
!..   for acceleration possibly use the inline-code of E3               
!..                                                                     
            A (L) = E3 (X (K) - X (I), Y (K) - Y (I), Z (K) - Z (I),    &
            M)                                                          
   10    END DO 
         L = L + 1 
!..                                                                     
!..   Set diagonal of G equal to zero                                   
!..                                                                     
         A (L) = 0.0D0 
   20 END DO 
      RETURN 
      END SUBROUTINE GAMMA3                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION APPRX3 (X0, Y0, Z0, NX, M, X, Y, Z, C) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Evaluation FUNCTION for the interpolation                      *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X0,Y0,Z0 :  location where function is to be evaluated         *      
! NX       :  number of nodes                                    *      
! M        :  derivative order                                   *      
! X,Y,Z    :  NX-vectors ..(1:NX); the coordinates of the nodes  *      
! C        :  vector of coefficients C(1:(NX + M*(M+1)*(M+2)/6)) *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! APPRX3   :  Approximate value at (X0,Y0,Z0)                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines:  none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (NX), Y (NX), Z (NX), C (NX + (M * (M + 1) * (M + 2)  &
      / 6) )                                                            
!..                                                                     
!..   for different M there are several cases:                          
!..   1. M = 1, 2, 3 ; especially coded, fast                           
!..   2. M > 3       ; each monomial is represented in the form         
!..                    (X**IX)*(Y**IY)*(Z**IZ), the evaluation is slow a
!..                    rounding error prone.                            
!..                                                                     
!..   the first polynomial is always  1                                 
!..                                                                     
      AP = C (NX + 1) 
      IF (M.EQ.1) GOTO 40 
      IF (M.EQ.2) THEN 
!..                                                                     
!..   remaining monomials of degree  1                                  
!..                                                                     
         AP = AP + C (NX + 2) * X0 + C (NX + 3) * Y0 + C (NX + 4)       &
         * Z0                                                           
      ELSEIF (M.EQ.3) THEN 
!..                                                                     
!..   remaining monomials of degree  2                                  
!..                                                                     
         AP = AP + (C (NX + 2) + C (NX + 5) * X0 + C (NX + 6) * Y0)     &
         * X0 + (C (NX + 3) + C (NX + 8) * Y0 + C (NX + 9) * Z0)        &
         * Y0 + (C (NX + 4) + C (NX + 7) * X0 + C (NX + 10) * Z0)       &
         * Z0                                                           
      ELSE 
!..                                                                     
!..   remaining monomials of degree <= M-1                              
!..                                                                     
         L = 1 
         DO 30 I = 1, M - 1 
            DO 20 IX = I, 0, - 1 
               DO 10 IY = I - IX, 0, - 1 
                  IZ = I - IX - IY 
                  L = L + 1 
                  IF (IX.NE.0.AND.IY.NE.0.AND.IZ.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) * (Y0**IY) *       &
                     (Z0**IZ)                                           
                  ELSEIF (IX.NE.0.AND.IY.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) * (Y0**IY) 
                  ELSEIF (IX.NE.0.AND.IZ.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) * (Z0**IZ) 
                  ELSEIF (IY.NE.0.AND.IZ.NE.0) THEN 
                     AP = AP + C (NX + L) * (Y0**IY) * (Z0**IZ) 
                  ELSEIF (IX.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) 
                  ELSEIF (IY.NE.0) THEN 
                     AP = AP + C (NX + L) * (Y0**IY) 
                  ELSE 
                     AP = AP + C (NX + L) * (Z0**IZ) 
                  ENDIF 
   10          END DO 
   20       END DO 
   30    END DO 
      ENDIF 
   40 CONTINUE 
!..                                                                     
!..   the G part of the system matrix                                   
!..                                                                     
!..   one might use the function E3(X,Y,Z) here, but this               
!..   would slow down the evaluation.                                   
!..   Hence this part is coded directly.                                
!..                                                                     
      DO 50 I = 1, NX 
         R = (X (I) - X0) **2 + (Y (I) - Y0) **2 + (Z (I) - Z0) **2 
         AP = AP + C (I) * (DSQRT (R) ** (2 * M - 3) ) 
   50 END DO 
      APPRX3 = AP 
      RETURN 
      END FUNCTION APPRX3                           
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION E3 (X, Y, Z, M) 
!                                                                       
!*****************************************************************      
!                                                                *      
! We evaluate the function F at (X,Y,Z).                         *      
!                                                                *      
! REMARK: the normalizing factors mentioned by MEINGUET [MEING79]*      
!         are not used as they turn out to be insignificant in   *      
!         the practical computations.                            *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X,Y,Z :  location where the evaluation takes place             *      
!                                                                *      
! M     :  derivative order                                      *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! E3    :  Functional value                                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines:  none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   Compute the kernel function                                       
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      R = DSQRT (X * X + Y * Y + Z * Z) 
      E3 = R** (2 * M - 3) 
      RETURN 
      END FUNCTION E3                               

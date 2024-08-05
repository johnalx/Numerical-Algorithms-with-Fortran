      SUBROUTINE CLINE (L, R, N, IA, COND, X, Y, Z, ZSUM, NA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!   This subroutine computes an estimate for the condition number*      
!   of a matrix A whose  LR decoposition is known, where L is a  *      
!   unit lower triangular and R is a nonsingular upper triangular*      
!   matrix.                                                      *      
!                                                                *      
!                                                                *      
!   INPUT PARAMETERS:                                            *      
!   =================                                            *      
!   L    : array L(1:IA,1:N) containing the unit diagonal lower  *      
!          triangular matrix of the LR decomposition of A        *      
!   R    : array R(1:IA,1:N) containing the nonsingular upper    *      
!          triangular matrix of the LR decomposition of A        *      
!   N    : order of the matrices L and R                         *      
!   IA   : leading dimension of the matrices L and R, as         *      
!          stipulated in the calling program                     *      
!                                                                *      
!                                                                *      
!   OUTPUT PARAMETER:                                            *      
!   =================                                            *      
!   COND : Estimate for the condition number of A                *      
!                                                                *      
!                                                                *      
!   AUXILIARY PARAMETERS:                                        *      
!   =====================                                        *      
!   X    : ] N-vectors ..(1:N)                                   *      
!   Y    : ]                                                     *      
!   Z    : ]                                                     *      
!   ZSUM : (N+1)-vector ZSUM(0:N)                                *      
!   NA   : array NA(1:IA,1:N)                                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!   Required subroutines: BACK, MXNORM, ZNORM                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!   Author    : Michaela Kisters                                 *      
!   Date      : 12.09.1990                                       *      
!   Source    : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION NA (IA, N), L (IA, N), R (IA, N), COND, X (N),    &
      Y (N), Z (N), ZSUM (0:N), MXNORM, ZNORM, KUNEND, SMI, SPL, V      
!                                                                       
!   For R transpose (=TRANS(R)) we determine X=X(I) with X(I)=+ or -1 an
!   Y=Y(I)=INV(TRANS(R))*X for I=1, ..., N, so that the 1-norm of Y beco
!   as large as possible.                                               
!                                                                       
      X (1) = 1.0D0 
      Y (1) = 1.0D0 / R (1, 1) 
      DO 10 I = 2, N 
         Y (I) = - R (1, I) * Y (1) / R (I, I) 
   10 END DO 
      DO 20 K = 2, N 
         V = 1.0D0 / R (K, K) 
         X (K) = Y (K) - V 
         Y (K) = Y (K) + V 
         SMI = DABS (X (K) ) 
         SPL = DABS (Y (K) ) 
         DO 30 I = K + 1, N 
            V = R (K, I) / R (I, I) 
            X (I) = Y (I) - V * X (K) 
            Y (I) = Y (I) - V * Y (K) 
            SMI = SMI + DABS (X (I) ) 
            SPL = SPL + DABS (Y (I) ) 
   30    END DO 
         IF (SMI.GT.SPL) THEN 
            DO 40 I = K, N 
               Y (I) = X (I) 
   40       END DO 
            X (K) = - 1.0D0 
         ELSE 
            X (K) = 1.0D0 
         ENDIF 
   20 END DO 
!                                                                       
!   Use backsubstitution to find Z with                                 
!   TRANS(L)*Z=Y.                                                       
!                                                                       
      CALL BACK (L, Y, N, IA, Z) 
!                                                                       
!   Estimate  KUNEND, the row sum norm of A                             
!                                                                       
      KUNEND = MXNORM (Z, N) / MXNORM (X, N) 
!                                                                       
!   Estimate COND(A)                                                    
!                                                                       
      COND = ZNORM (L, R, N, IA, ZSUM, NA) * KUNEND 
!                                                                       
      RETURN 
      END SUBROUTINE CLINE                          
!                                                                       
!                                                                       
      SUBROUTINE BACK (L, B, N, IL, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!   Solving a triangular system TRANS(L)*X=B by backsubstitution *      
!                                                                *      
!   INPUT PARAMETERS:                                            *      
!   =================                                            *      
!   L    : array L(1:IL,1:N) contaning the entries of the unit   *      
!          diagonal lower tringular matrix L                     *      
!   B    : N-vector B(1:N), the right hand side of the system    *      
!   N    : order of the matrix L                                 *      
!   IL   : leading dimension of L, as stipulated in the calling  *      
!          program                                               *      
!                                                                *      
!                                                                *      
!   OUTPUT PARAMETER:                                            *      
!   =================                                            *      
!   X    : N-vector X(1:N), the solution vector                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!   Required subroutines: none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!   Author    : Michaela Kisters                                 *      
!   Date      : 09.12.1990                                       *      
!   Sourcee   : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION L (IL, N), B (N), X (N), SUM 
!                                                                       
      X (N) = B (N) / L (N, N) 
      DO 10 K = N - 1, 1, - 1 
         SUM = 0.0D0 
         DO 20 J = K + 1, N 
            SUM = SUM + L (J, K) * X (J) 
   20    END DO 
         X (K) = 1.0D0 / L (K, K) * (B (K) - SUM) 
   10 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE BACK                           
!                                                                       
!                                                                       
      FUNCTION MXNORM (X, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!   Calculate the maximum norm MXNORM of a vector X              *      
!                                                                *      
!   INPUT PARAMETERS:                                            *      
!   =================                                            *      
!   X    : N-vector X(1:N)                                       *      
!   N    : dimension of X                                        *      
!                                                                *      
!                                                                *      
!   OUTPUT PARAMETER:                                            *      
!   =================                                            *      
!   MXNORM : Maximum norm of X                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!   Required subroutines: none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!   Author    : Michaela Kisters                                 *      
!   Date      : 12.09.1990                                       *      
!   Source    : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X (N), MXNORM 
!                                                                       
      MXNORM = DABS (X (1) ) 
      DO 10 I = 1, N 
         IF (DABS (X (I) ) .GT.MXNORM) MXNORM = DABS (X (I) ) 
   10 END DO 
!                                                                       
      RETURN 
      END FUNCTION MXNORM                           
!                                                                       
!                                                                       
      FUNCTION ZNORM (L, R, N, IA, ZSUM, NA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!   Computes the row sum norm of a matrix A which is given by its*      
!   LR decomposition by a unit diagonal lower triangular matrix L*      
!   and a nonsingular upper tringular matrix R.                  *      
!                                                                *      
!                                                                *      
!   INPUT PARAMETERS:                                            *      
!   =================                                            *      
!   L    : array L(1:IA,1:N) containing the unit diagonal lower  *      
!          triangular matrix of the LR decomposition of A        *      
!   R    : array R(1:IA,1:N) containing the nonsingular upper    *      
!          triangular matrix of the LR decomposition of A        *      
!   N    : order of the matrices L and R                         *      
!   IA   : leading dimension of the matrices L and R, as         *      
!          stipulated in the calling program                     *      
!                                                                *      
!                                                                *      
!   OUTPUT PARAMETER:                                            *      
!   =================                                            *      
!   ZNORM : Row sum norm of the matrix A                         *      
!                                                                *      
!                                                                *      
!   AUXILIARY PARAMETERS:                                        *      
!   =====================                                        *      
!   ZSUM : (N+1)-vector ZSUM(0:N)                                *      
!   NA   : array NA(1:IA,1:N)                                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!   Required subroutines: none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!   Author    : Michaela Kisters                                 *      
!   Date      : 12.09.1990                                       *      
!   Source    : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION ZSUM (0:N), ZNORM, L (IA, N), R (IA, N), NA (IA,  &
      N)                                                                
!                                                                       
      ZSUM (0) = - 99.0D0 
      DO 10 I = 1, N 
         ZSUM (I) = 0.0D0 
         DO 20 J = 1, N 
            NA (I, J) = 0.0D0 
            DO 30 K = 1, I 
               NA (I, J) = NA (I, J) + L (I, K) * R (K, J) 
   30       END DO 
            ZSUM (I) = ZSUM (I) + DABS (NA (I, J) ) 
   20    END DO 
         IF (ZSUM (I) .GT.ZSUM (I - 1) ) ZNORM = ZSUM (I) 
   10 END DO 
!                                                                       
      RETURN 
      END FUNCTION ZNORM                            

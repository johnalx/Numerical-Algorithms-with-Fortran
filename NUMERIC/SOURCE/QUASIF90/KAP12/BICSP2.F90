      SUBROUTINE BICSP2 (N, M, A, X, Y, F, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Determination of bicubic splines.                             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    : number of X-intervals                                  *      
!  M    : number of Y-intervals                                  *      
!  A    : 4-dimensional array A(0:N,0:M,0:KDIM,0:LDIM); contains *      
!         the spline coefficients. On call, A(I,J,0,0) must con- *      
!         tain the functional values U(I,J).                     *      
!         BICSP1 determines all other A(I,J,K,L) for I=0 to N-1  *      
!         and J=0 to M-1. Elements A(N,M,K,L), that are not      *      
!         assigned any values on call, remain unassigned.        *      
!  X    : (N+1)-vector X(0:N) contains the endpoints of the      *      
!         X-intervals                                            *      
!  Y    : (N+1)-vector Y(0:M) contains the endpoints of the      *      
!         Y-intervals                                            *      
!  F    : auxiliary vector F(1:9*MAX(N,M)-8)                     *      
!  IERR : is initially set to 0. Will be changed to be nonzero,  *      
!         if the algorithm detects an error. If errors occur,    *      
!         the program does not complete the computations.        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BIC2S1, BIC2S2, BIC2S3, BICSP1          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), Y (0:M), F ( * ) 
!                                                                       
!*  step 1                                                              
!                                                                       
      CALL BIC2S1 (N, M, A, X, F, JERR) 
      IERR = JERR 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 2                                                              
!                                                                       
      CALL BIC2S2 (N, M, A, Y, F, JERR) 
      IERR = JERR + 1 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 3                                                              
!                                                                       
      CALL BIC2S3 (N, M, A, X, F, JERR) 
      IERR = JERR + 4 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 4 to 12                                                        
!                                                                       
      CALL BICSP1 (N, M, A, X, Y, F, JERR) 
      IERR = JERR + 4 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  algorithm has run successfully                                      
!                                                                       
      IERR = 0 
      RETURN 
      END SUBROUTINE BICSP2                         
!                                                                       
!                                                                       
      SUBROUTINE BIC2S1 (N, M, A, X, H, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  step 1:                                                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), H (0:N - 1) 
!                                                                       
!*  determine H(0),H(1),H(N-2),H(N-1), test strict monotonicity         
!                                                                       
      DO 102 K = 0, 1 
         DO 101 L = 0, N - 2, N - 2 
            I = K + L 
            H (I) = X (I + 1) - X (I) 
            IF (H (I) .LE.0.0D0) GOTO 900 
  101    END DO 
  102 END DO 
!                                                                       
!*  determine boundary values and store them in A                       
!                                                                       
      DO 103 J = 0, M 
!                                                                       
!*  converted formulas                                                  
!                                                                       
         A (0, J, 1, 0) = (A (1, J, 0, 0) - A (0, J, 0, 0) ) * (1.0D0 / &
         H (0) + 0.5D0 / (H (0) + H (1) ) ) - (A (2, J, 0, 0) - A (1, J,&
         0, 0) ) * H (0) / (H (1) * 2.0D0 * (H (0) + H (1) ) )          
         A (N, J, 1, 0) = (A (N, J, 0, 0) - A (N - 1, J, 0, 0) )        &
         * (0.5D0 / (H (N - 2) + H (N - 1) ) + 1.0D0 / H (N - 1) )      &
         - (A (N - 1, J, 0, 0) - A (N - 2, J, 0, 0) ) * H (N - 1)       &
         / (2.0D0 * (H (N - 2) + H (N - 1) ) * H (N - 2) )              
!                                                                       
  103 END DO 
      IERR = 0 
      RETURN 
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC2S1                         
!                                                                       
!                                                                       
      SUBROUTINE BIC2S2 (N, M, A, Y, H, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  step 2:                                                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), Y (0:M), H (0:M - 1) 
!                                                                       
!*  determine H(0),H(1),H(M-2),H(M-1), test strict monotonicity         
!                                                                       
      DO 102 K = 0, 1 
         DO 101 L = 0, M - 2, M - 2 
            J = K + L 
            H (J) = Y (J + 1) - Y (J) 
            IF (H (J) .LE.0.0D0) GOTO 900 
  101    END DO 
  102 END DO 
!                                                                       
!*  determine boundary values and store them in A                       
!                                                                       
      DO 103 I = 0, N 
!                                                                       
!*  converted formulas                                                  
!                                                                       
         A (I, 0, 0, 1) = (A (I, 1, 0, 0) - A (I, 0, 0, 0) ) * (1.0D0 / &
         H (0) + 0.5D0 / (H (0) + H (1) ) ) - (A (I, 2, 0, 0) - A (I, 1,&
         0, 0) ) * H (0) / (H (1) * 2.0D0 * (H (0) + H (1) ) )          
         A (I, M, 0, 1) = (A (I, M, 0, 0) - A (I, M - 1, 0, 0) )        &
         * (0.5D0 / (H (M - 2) + H (M - 1) ) + 1.0D0 / H (M - 1) )      &
         - (A (I, M - 1, 0, 0) - A (I, M - 2, 0, 0) ) * H (M - 1)       &
         / (2.0D0 * (H (M - 2) + H (M - 1) ) * H (M - 2) )              
!                                                                       
  103 END DO 
      IERR = 0 
      RETURN 
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC2S2                         
!                                                                       
!                                                                       
      SUBROUTINE BIC2S3 (N, M, A, X, H, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  step 3:                                                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), H (0:N - 1) 
!                                                                       
!*  determine H(0),H(1),H(N-2),H(N-1), test strict monotonicity         
!                                                                       
      DO 102 K = 0, 1 
         DO 101 L = 0, N - 2, N - 2 
            I = K + L 
            H (I) = X (I + 1) - X (I) 
            IF (H (I) .LE.0.0D0) GOTO 900 
  101    END DO 
  102 END DO 
!                                                                       
!*  determine boundary values and store them in A                       
!                                                                       
      DO 103 J = 0, M, M 
!                                                                       
!*  converted formulas                                                  
!                                                                       
         A (0, J, 1, 1) = (A (1, J, 0, 1) - A (0, J, 0, 1) ) * (1.0D0 / &
         H (0) + 0.5D0 / (H (0) + H (1) ) ) - (A (2, J, 0, 1) - A (1, J,&
         0, 1) ) * H (0) / (H (1) * 2.0D0 * (H (0) + H (1) ) )          
         A (N, J, 1, 1) = (A (N, J, 0, 1) - A (N - 1, J, 0, 1) )        &
         * (0.5D0 / (H (N - 2) + H (N - 1) ) + 1.0D0 / H (N - 1) )      &
         - (A (N - 1, J, 0, 1) - A (N - 2, J, 0, 1) ) * H (N - 1)       &
         / (2.0D0 * (H (N - 2) + H (N - 1) ) * H (N - 2) )              
!                                                                       
  103 END DO 
      IERR = 0 
      RETURN 
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC2S3                         

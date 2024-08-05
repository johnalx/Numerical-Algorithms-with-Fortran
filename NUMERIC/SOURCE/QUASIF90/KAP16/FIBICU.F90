![  {Riemann Double Integrals using Bicubic Splines}                    
![  {Riemann Double Integrals using Bicubic Splines}*)                  
      SUBROUTINE FIBICU (N, M, A, X, Y, VALUE) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine determines a double integral of a spline      *      
!  function over the rectangle [X(0), X(N)] x [Y(0), Y(M)],      *      
!  which is the complete domain of the spline.                   *      
!  This program is faster than FIBIC2.                           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: FIBIC1                                  *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), Y (0:M) 
!                                                                       
      S = 0.0D0 
      DO 102 I = 0, N - 1 
         DO 101 J = 0, M - 1 
            S = S + FIBIC1 (N, M, A, I, J, X (I + 1) - X (I), Y (J + 1) &
            - Y (J) )                                                   
  101    END DO 
  102 END DO 
      VALUE = S 
      RETURN 
      END SUBROUTINE FIBICU                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION FIBIC1 (N, M, A, I, J, XI, ETA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The function FIBIC1 determines a double integral of a spline  *      
!  function over the rectangle [0, XI] x [0, ETA]. XI and ETA    *      
!  are relative coordinates in the rectangle I, J.               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  autor     : Eberhard Heyne                                    *      
!  datum     : 02.15.1983                                        *      
!  sourcecode: FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM) 
      DIMENSION XIP (0:4), ETAP (0:4) 
      DATA XIP (0) / 1.0D0 / 
      DATA ETAP (0) / 1.0D0 / 
      S = 0.0D0 
      DO 101 K = 1, 4 
         XIP (K) = XIP (K - 1) * XI 
         ETAP (K) = ETAP (K - 1) * ETA 
  101 END DO 
      DO 103 K = 0, 3 
         DO 102 L = 0, 3 
            S = S + A (I, J, K, L) * XIP (K + 1) * ETAP (L + 1) / DBLE (&
            (K + 1) * (L + 1) )                                         
  102    END DO 
  103 END DO 
      FIBIC1 = S 
      RETURN 
      END FUNCTION FIBIC1                           

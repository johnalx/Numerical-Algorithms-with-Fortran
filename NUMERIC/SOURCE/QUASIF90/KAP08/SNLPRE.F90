      SUBROUTINE SNLPRE (X, W, PHI, DVT, JNDVT, C, LDA, M, N, EPS, F, A) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The SUBROUTINE SNLPRE computes the Jacobi matrix (as it is    *      
!  needed in SUBROUTINE SNLFIT), either using a subroutine,      *      
!  denoted by DVT in the program that is provided by the user    *      
!  and which computes the partial derivatives, or by central     *      
!  difference quotients.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X, PHI, DVT, JNDVT, C, LDA, M, N      same as SNLFIT          *      
!                                                                *      
!  W     a vector which contains the square roots of the weights *      
!                                                                *      
!  EPS   the machine constant                                    *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!                                                                *      
!  F     (N+1)-vector F(0:N), which is needed when calling the   *      
!         user supplied SUBROUTINE DVT                           *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!                                                                *      
!  A     2-dim. array A(0:LDA,0:N) containing the Jacobi matrix  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none, except possibly DVT, which must  *      
!                         be user supplied if JVDT = 0           *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Ilona Westermann                                   *      
!  date     : 09.01.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER JNDVT, LDA, M, N 
      DIMENSION X (0:M), W (0:M), C (0:N), F (0:N), A (0:LDA, 0:N) 
!                                                                       
!  determine the partial derivatives using SUBROUTINE DVT               
!                                                                       
      IF (JNDVT.EQ.0) THEN 
         DO 10 I = 0, M 
            CALL DVT (X (I), C, N, F) 
            DO 10 K = 0, N 
               A (I, K) = F (K) 
   10    CONTINUE 
!                                                                       
!  approximate the partial derivatives by central                       
!  difference quotients                                                 
!                                                                       
      ELSE 
         FACTOR = EPS** (1.0D0 / 3.0D0) 
         DO 20 K = 0, N 
            IF (C (K) .EQ.0.0D0) THEN 
               HK = FACTOR 
            ELSE 
               HK = FACTOR * DABS (C (K) ) 
            ENDIF 
            ZHK = 1.0D0 / (2.0D0 * HK) 
            DO 20 I = 0, M 
               C (K) = C (K) + HK 
               DIFQUO = PHI (C, N, X (I) ) 
               C (K) = C (K) - 2.0D0 * HK 
               A (I, K) = (DIFQUO - PHI (C, N, X (I) ) ) * ZHK 
               C (K) = C (K) + HK 
   20    CONTINUE 
      ENDIF 
      DO 30 I = 0, M 
         DO 30 K = 0, N 
            A (I, K) = A (I, K) * W (I) 
   30 CONTINUE 
      RETURN 
      END SUBROUTINE SNLPRE                         

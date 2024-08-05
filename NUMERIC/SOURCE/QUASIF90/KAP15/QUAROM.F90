      SUBROUTINE QUAROM (A, B, EPS, N, H, FCT, EL, RESULT, ERREST, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This subroutine determines an approximation for the        *      
!     integral of the FUNCTION FCT(X) over the interval [A,B]    *      
!     using the ROMBERG-method.                                  *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     A,B    - the interval endpoints.                           *      
!     EPS    - accuracy bound for the error estimate.            *      
!     N      - maximum number of rows and columns of the         *      
!              ROMBERG scheme. (N > 1)                           *      
!     H      - starting step size for which the following must   *      
!              hold:                                             *      
!                    H = (B-A) / K  for  K a positive integer.   *      
!              If H was chosen wrongly, H is interally set to    *      
!              equal (B-A) (without any specific error message). *      
!     FCT    - function to be integrated.                        *      
!              It has to be provided by the user in the following*      
!              format:                                           *      
!                     DOUBLE PRECISION FUNCTION  FCT (x).        *      
!              The function has to be defined as EXTERNAL in the *      
!              calling program.                                  *      
!     EL     - auxiliary vector of length N at least.            *      
!              In EL the current row of the ROMBERG scheme       *      
!              is stored.                                        *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     RESULT - approximate value for the integral at the end of  *      
!              the procedure                                     *      
!     ERREST - error estimate for the approximate value RESULT   *      
!     N      - number of rows (columns) of the ROMBERG scheme    *      
!              that were actually determined                     *      
!     H      - step size at end of calculations                  *      
!     IERR   - error parameter,                                  *      
!                IERR = 0 : everything o.k. ERREST < EPS         *      
!                IERR = 1 : incorrect input parameters :         *      
!                            N < 1 or EPS < 0.0  .               *      
!                IERR = 2 : required accuracy was not achieved   *      
!                           after N steps, i.e., ERREST > EPS .  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Richard Reuter  (1983, FORTRAN IV)                 *      
!  editor   : Gisela Engeln-Muellges (1988)                      *      
!  editor   : Norbert Vogt                                       *      
!  date     : 01.31.1990                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  declarations                                                         
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION EL (N) 
!                                                                       
!  testing the input data                                               
!                                                                       
      IERR = 1 
      IF (N.LE.1.OR.EPS.LE.0.0D0) RETURN 
      IF (A.EQ.B) THEN 
         RESULT = 0.0D0 
         N = 0 
         ERREST = 0.0D0 
         IERR = 0 
         H = 0.0D0 
         RETURN 
      ENDIF 
!                                                                       
!  determine the starting step size;                                    
!  determine the number N0 of sub-intervals                             
!                                                                       
      H = DMIN1 (DABS (H), DABS (B - A) ) 
      IF (B.LT.A) H = - H 
      IF (H.EQ.0.0D0) H = B - A 
      N0 = (B - A + 0.5D0 * H) / H 
!                                                                       
!  determine the first row of the ROMBERG scheme                        
!                                                                       
      IERR = 0 
      EL (1) = 0.5D0 * (FCT (A) + FCT (B) ) 
      IF (N0.NE.1) THEN 
         DO 10 L = 2, N0 
            EL (1) = EL (1) + FCT (A + (L - 1) * H) 
   10    END DO 
      ENDIF 
!                                                                       
!  approximate integral from the first ROMBERG scheme row               
!                                                                       
      EL (1) = H * EL (1) 
!                                                                       
!  determine further rows of the ROMBERG scheme                         
!                                                                       
      DO 20 K = 2, N 
         EL (K) = 0.0D0 
         H = H * 0.5D0 
         EL1 = EL (1) 
         EL (1) = 0.0D0 
         DO 30 L = 1, N0 
            EL (1) = EL (1) + FCT (A + (2 * L - 1) * H) 
   30    END DO 
         EL (1) = EL (1) * H + EL1 * 0.5D0 
         N0 = 2 * N0 
!                                                                       
!  determine the linear combination in the K-th row                     
!                                                                       
         MM = 1 
         DO 40 M = 2, K 
            MM = MM * 4 
            EL2 = EL (M) 
            EL (M) = (MM * EL (M - 1) - EL1) / (MM - 1) 
            EL1 = EL2 
   40    END DO 
!                                                                       
!  determine an estimate for the error of RESULT                        
!                                                                       
         ERREST = DABS (EL (K) - EL (K - 1) ) 
!                                                                       
!  decide upon stopping the procedure                                   
!                                                                       
         IF (ERREST.LT.EPS) THEN 
            RESULT = EL (K) 
            N = K 
            RETURN 
         ENDIF 
   20 END DO 
      IERR = 2 
      RESULT = EL (N) 
      RETURN 
      END SUBROUTINE QUAROM                         

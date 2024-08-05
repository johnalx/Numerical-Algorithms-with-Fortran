      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Testprogram for subroutines GAUSS, GAUSSP, GAUSSS, HACOND,    *      
!  POSTIT.                                                       *      
!                                                                *      
!  Solve linear system   A * X = Y  with  Hilbert matrix         *      
!                                                                *      
!  A(I,K)=1./(I+K-1) ,  I,K = 1, .., N  for  N = 1, ..., 10      *      
!                                                                *      
!  and  right hand side Y(I)=A(I,1)+A(I,2)+...+A(I,N) , I=1,..,N.*      
!                                                                *      
!  Solution:  X(I) = 1. , I = 1, ..., N.                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (1:10, 1:10), A0 (1:10, 1:10), Y (1:10),        &
      X (1:10), D (1:10), Z (1:10), RS (1:10), R (1:10)                 
      INTEGER IPIVOT (1:10) 
      CHARACTER FORM * 30, TEX1 * 2 
      EPS = 0.1D-12 
      DO 10 N = 1, 10 
         DO 20 I = 1, N 
            SUM = 0.D0 
            DO 30 K = 1, N 
               A (I, K) = 1.D0 / (I + K - 1.D0) 
               A0 (I, K) = A (I, K) 
               SUM = SUM + A (I, K) 
   30       END DO 
            Y (I) = SUM 
   20    END DO 
         WRITE (TEX1, '(I2)') N 
         FORM = '(1X,'//TEX1//'(F5.3,1X),2X,F5.3)' 
         DO 40 I = 1, N 
            WRITE ( *, FORM) (A (I, K), K = 1, N), Y (I) 
   40    END DO 
         CALL GAUSS (N, A, 10, Y, X, MARKE, D, IPIVOT) 
         WRITE ( * , '(//,1X,''SOLUTION:'')') 
         WRITE ( * , * ) '*********' 
         WRITE ( *, * ) 
         WRITE ( * , '(/,1X,''MARK = '',I3,//)') MARKE 
         IF (MARKE.NE.0) THEN 
            DO 50 I = 1, N 
               WRITE ( * , '(1X,D20.14,10X,D20.14)') X (I) , DABS (X (I)&
               - 1.)                                                    
   50       END DO 
            CALL HACOND (N, A0, A, 10, MARKE, HKOND) 
            WRITE ( * , '(//,1X,''CONDITION NUMBER: '',D20.14)') HKOND 
            CALL POSTIT (N, A0, A, 10, IPIVOT, Y, X, EPS, 10, ITANZ,    &
            IFEHL, Z, R, RS)                                            
            WRITE ( * , '(//,1X,''BREAK-OFF CONDITION: '',I3)') IFEHL 
            IF (IFEHL.EQ.0) THEN 
               WRITE ( *, 900) ITANZ, EPS 
               DO 60 I = 1, N 
                  WRITE ( * , '(1X,D20.14,10X,D20.14)') X (I) , DABS (X &
                  (I) - 1.)                                             
   60          END DO 
            ENDIF 
         ENDIF 
         WRITE ( * , '(1H1)') 
   10 END DO 
      STOP 
  900 FORMAT(//,1X,'NUMBER OF ITERATIONS: ',I3,10X,'FAULIY BOUND: ',    &
     &       D20.14,//)                                                 
      END PROGRAM TEST                              

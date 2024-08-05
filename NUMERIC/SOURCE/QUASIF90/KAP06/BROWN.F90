![           {Brown's Method for Nonlinear Systems}*)                   
      SUBROUTINE BROWN (N, FCT, X0, EPS, EPSM, IOUT, MAXIT, IHF, LDIHF, &
      HF, LDHF, X1, NUMIT, IERR)                                        
!                                                                       
!*****************************************************************      
!                                                                *      
! SUBROUTINE BROWN determines the zero of a nonlinear system of  *      
! equations consisting of N equations in N unknowns by applying  *      
! the Brown-method.                                              *      
!                                                                *      
! In order to perfornm one iteration step with Brown's algorithm *      
! the SUBROUTINE ITER4 is called. The iteration is continued     *      
! until the given maximal number of iteration steps is reached.  *      
! The iteration is stopped if one of the following break-off     *      
! criteria is met:                                               *      
! - if the relative change of two successive iterative solutions *      
!   is smaller than EPS.                                         *      
! - if the functional value is smaller or equal to EPSM.         *      
! - if the limiting accuracy is reached.                         *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! N       : number of equations and unknowns                     *      
! FCT     : function subroutine, given as follows:               *      
!                 SUBROUTINE FCT (K,X,F)                         *      
!                 DOUBLE PRECISION X(2)                          *      
!                 GO TO (1,2) K                                  *      
!               1 F = F1(X(1),X(2))                              *      
!                 RETURN                                         *      
!               2 F = F2(X(1),X(2))                              *      
!                 RETURN                                         *      
!                 END                                            *      
! X0      : N-vector X0(1:N); starting vector with N components  *      
! EPS     : desired accuracy                                     *      
! EPSM    : machine constant                                     *      
! IOUT    : IOUT=1 after each iteration step the program puts    *      
!                  out the difference of the last approximate    *      
!                  solutions, the current approximate solution   *      
!                  and the functional value                      *      
!           IOUT=0 no output is provided, except upon            *      
!                  convergence                                   *      
! MAXIT   : maximum number of iteration steps allowed            *      
! IHF     : 2-dim. array IHF(1:LDIHF,1:N+1); ) auxiliary arrays  *      
! HF      : 2-dim. array   HF(1:LDHF,1:N+3); ) used as           *      
!                                            ) intermediate      *      
!                                            ) storage space     *      
! LDIHF   : leading dimension of IHF as defined in the calling   *      
!           program. LDIHF >= N                                  *      
! LDHF    : leading dimension of HF as defined in the calling    *      
!           program. LDHF >= N                                   *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X1      : N-vector X1(1:N); the approximate solution vector    *      
! NUMIT   : number of iterations executed                        *      
! IERR    : error parameter                                      *      
!           =0; no error, zero has been found                    *      
!           =1; after MAXIT steps the desired accuracy was not   *      
!               reached                                          *      
!           =2; singular matrix                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: ITER4, SUBST                            *      
!                                                                *      
!                                                                *      
!  sources : Brown K.M.: A quadratically convergent Newton-like  *      
!            method based upon Gaussian elimination,             *      
!            Siam J.Numer.Anal.vol 6 (1969), 560 - 569           *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Johannes Karfusehr                                 *      
!  Editor   : Thomas Eul                                         *      
!  Date     : 08.21.1985                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!   declarations                                                        
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      EXTERNAL FCT 
      INTEGER N, IERR, LDIHF, LDHF 
      DIMENSION X0 (N), X1 (N) 
      DIMENSION IHF (LDIHF, N + 1), HF (LDHF, N + 3) 
!                                                                       
!   local variables                                                     
!                                                                       
      LOGICAL SING 
!                                                                       
!   initializing                                                        
!                                                                       
      SING = .FALSE. 
      ICRIT = 0 
      DELTA0 = 1.0D-2 
      DO 10 J = 1, N 
         HF (J, N + 2) = X0 (J) 
   10 END DO 
      IF (IOUT.EQ.1) WRITE ( *, 1000) 
!                                                                       
!   iteration                                                           
!                                                                       
      DO 20 M = 1, MAXIT 
         CALL ITER4 (N, FCT, EPSM, IHF, LDIHF, HF, LDHF, X1, SING) 
         NUMIT = M 
!                                                                       
!   after each iteration step we output the current                     
!   starting vector and the next iterate                                
!                                                                       
         IF (IOUT.EQ.1) THEN 
            WRITE ( *, 3000) M 
            IF (.NOT.SING) THEN 
!                                                                       
               WRITE ( *, 4000) 
               DO 30 J = 1, N 
                  CALL FCT (J, X1, FVAL) 
                  WRITE ( *, 5000) X1 (J) - HF (J, N + 2), J, X1 (J),   &
                  FVAL                                                  
   30          END DO 
            ELSE 
               WRITE ( *, 7000) 
            ENDIF 
         ENDIF 
!                                                                       
         IF (.NOT.SING) THEN 
!                                                                       
!   test the break-off criterion                                        
!                                                                       
!   test the relative change of the iterates                            
!                                                                       
            DO 40 I = 1, N 
               RELF = (X1 (I) - HF (I, N + 2) ) / (HF (I, N + 2)        &
               + EPS)                                                   
               IF (DABS (RELF) .GE.EPS) THEN 
                  GOTO 41 
               ENDIF 
   40       END DO 
            ICRIT = 1 
            GOTO 21 
   41       CONTINUE 
!                                                                       
!   test the functional value                                           
!                                                                       
            DO 50 I = 1, N 
               CALL FCT (I, X1, FVAL) 
               IF (DABS (FVAL) .GT.EPSM) GOTO 51 
   50       END DO 
            ICRIT = 2 
            GOTO 21 
   51       CONTINUE 
!                                                                       
!   test the limiting accuracy                                          
!                                                                       
            DELTA1 = DABS (X1 (1) - HF (1, N + 2) ) 
            DO 60 I = 2, N 
               DELTA1 = DMAX1 (DELTA1, DABS (X1 (I) - HF (I, N + 2) ) ) 
   60       END DO 
            IF (DELTA1.LE.1.0D-3) THEN 
               IF (DELTA0.LE.DELTA1) THEN 
                  ICRIT = 3 
                  GOTO 21 
               ENDIF 
            ENDIF 
            DELTA0 = DELTA1 
            IF (NUMIT.LT.MAXIT) THEN 
               DO 70 I = 1, N 
                  HF (I, N + 2) = X1 (I) 
   70          END DO 
            ENDIF 
         ELSE 
            GOTO 21 
         ENDIF 
   20 END DO 
   21 CONTINUE 
      IF (SING) THEN 
         IERR = 2 
      ELSE 
         IF (ICRIT.EQ.0) THEN 
            IERR = 1 
         ELSE 
            IERR = 0 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT ('1') 
 3000 FORMAT (1X,I4,'-TH ITERATION STEP :') 
 4000 FORMAT (1X,6X,'DIFFERENCE',7X,2X,'KOMP',2X,6X,'APPROXIMATION',7X, &
     &        2X,4X,'FUNCTIONAL VALUE')                                 
 5000 FORMAT (1X,D22.15,2X,I4,2X,D22.15,2X,D22.15) 
 7000 FORMAT(1X,'THE JACOBI MATRIX IS SINGULAR') 
!                                                                       
      RETURN 
      END SUBROUTINE BROWN                          
!                                                                       
!                                                                       
      SUBROUTINE ITER4 (N, FCT, EPSM, IHF, LDIHF, HF, LDHF, X1, SING) 
!                                                                       
!*****************************************************************      
!                                                                *      
! SUBROUTINE ITER4 computes an approximation by Brown's          *      
! algorithm.                                                     *      
!                                                                *      
!                                                                *      
! INPUT PARAMETER:                                               *      
! ================                                               *      
! N       : number of equations or unknowns                      *      
! FCT     : function subroutine (compare with SUBROUTINE BROWN)  *      
! EPSM    : machine constant                                     *      
! IHF     : 2-dim. array IHF(1:LDIHF,1:N+1); ) auxiliary array   *      
! HF      : 2-dim. array   HF(1:LDHF,1:N+3); ) providing storage *      
!                                            ) space             *      
! LDIHF   : leading dimension of IHF as defined in the calling   *      
!           program. LDIHF >= N                                  *      
! LDHF    : leading dimension of HF as defined in the calling    *      
!           program. LDHF >= N                                   *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X1      : N-vector X1(1:N); the approximation vector           *      
! SING    : error parameter; indicates whether the Jacobi matrix *      
!           is singular or not                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SUBST                                   *      
!                                                                *      
!                                                                *      
!  sources : Brown K.M.: A quadratically convergent Newton-like  *      
!            method based upon Gaussian elimination,             *      
!            Siam J.Numer.Anal.,vol 6 (1969), 560 - 569          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Johannes Karfusehr                                 *      
!  editor   : Thomas Eul                                         *      
!  date     : 08.21.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!     declarations                                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      LOGICAL SING 
      INTEGER N, LDIHF, LDHF 
      DIMENSION IHF (LDIHF, N + 1), HF (LDHF, N + 3), X1 (N) 
!                                                                       
!   initializing                                                        
!                                                                       
      DO 10 J = 1, N 
         IHF (1, J) = J 
         X1 (J) = HF (J, N + 2) 
   10 END DO 
!                                                                       
!   linearization of the K-th coordinate function                       
!                                                                       
      DO 20 K = 1, N 
         ICOUNT = 0 
         FACTOR = 0.001D0 
         DO 30 J = 1, 3 
            IF (K.GT.1) CALL SUBST (K, N, IHF, LDIHF, HF, LDHF, X1) 
            CALL FCT (K, X1, F) 
!                                                                       
!   determining the I-th discretization stepsize                        
!   and the I-th difference quotient                                    
!                                                                       
            DO 40 I = K, N 
               ITEMP = IHF (K, I) 
               HOLD = X1 (ITEMP) 
               H = FACTOR * HOLD 
               IF (DABS (H) .LE.EPSM) H = 0.001D0 
               X1 (ITEMP) = HOLD+H 
               IF (K.GT.1) CALL SUBST (K, N, IHF, LDIHF, HF, LDHF, X1) 
               CALL FCT (K, X1, FPLUS) 
               X1 (ITEMP) = HOLD 
               HF (ITEMP, N + 3) = (FPLUS - F) / H 
               IF (DABS (HF (ITEMP, N + 3) ) .LE.EPSM) THEN 
                  ICOUNT = ICOUNT + 1 
               ELSE 
                  IF (DABS (F / HF (ITEMP, N + 3) ) .GE.1.0D20) ICOUNT =&
                  ICOUNT + 1                                            
               ENDIF 
   40       END DO 
            IF (ICOUNT.LE.N - K) THEN 
               SING = .FALSE. 
               GOTO 31 
            ELSE 
               SING = .TRUE. 
               FACTOR = FACTOR * 10.0D0 
               ICOUNT = 0 
            ENDIF 
   30    END DO 
   31    CONTINUE 
!                                                                       
         IF (.NOT.SING) THEN 
            IF (K.LT.N) THEN 
               KMAX = IHF (K, K) 
!                                                                       
!   determining the difference quotient of largest magnitude            
!                                                                       
               DERMAX = DABS (HF (KMAX, N + 3) ) 
               KPLUS = K + 1 
               DO 50 I = KPLUS, N 
                  JSUB = IHF (K, I) 
                  TEST = DABS (HF (JSUB, N + 3) ) 
                  IF (TEST.LT.DERMAX) THEN 
                     IHF (KPLUS, I) = JSUB 
                  ELSE 
                     IHF (KPLUS, I) = KMAX 
                     KMAX = JSUB 
                  ENDIF 
   50          END DO 
               IF (DABS (HF (KMAX, N + 3) ) .LE.EPSM) SING = .TRUE. 
               IHF (K, N + 1) = KMAX 
               IF (.NOT.SING) THEN 
                  HF (K, N + 1) = 0.0D0 
!                                                                       
!   solving the K-th equation for XMAX                                  
!                                                                       
                  DO 60 J = KPLUS, N 
                     JSUB = IHF (KPLUS, J) 
                     HF (K, JSUB) = - HF (JSUB, N + 3) / HF (KMAX, N +  &
                     3)                                                 
                     HF (K, N + 1) = HF (K, N + 1) + HF (JSUB, N + 3)   &
                     * X1 (JSUB)                                        
   60             END DO 
                  HF (K, N + 1) = (HF (K, N + 1) - F) / HF (KMAX, N + 3)&
                  + X1 (KMAX)                                           
               ELSE 
                  GOTO 21 
               ENDIF 
            ELSE 
!                                                                       
!   solving the N-th coordinate function by use of the                  
!   discrete Newton-method for one variable                             
!                                                                       
               IF (DABS (HF (ITEMP, N + 3) ) .LE.EPSM) THEN 
                  SING = .TRUE. 
               ELSE 
                  HF (K, N + 1) = 0.0D0 
                  KMAX = ITEMP 
                  HF (K, N + 1) = (HF (K, N + 1) - F) / HF (KMAX, N + 3)&
                  + X1 (KMAX)                                           
               ENDIF 
            ENDIF 
         ELSE 
            GOTO 21 
         ENDIF 
   20 END DO 
   21 CONTINUE 
      IF (.NOT.SING) THEN 
!                                                                       
!   determining of the approximate solution by backsubstitution         
!                                                                       
         X1 (KMAX) = HF (N, N + 1) 
         IF (N.GT.1) CALL SUBST (N, N, IHF, LDIHF, HF, LDHF, X1) 
      ENDIF 
      RETURN 
      END SUBROUTINE ITER4                          
!                                                                       
!                                                                       
      SUBROUTINE SUBST (K, N, IHF, LDIHF, HF, LDHF, X1) 
!                                                                       
!*****************************************************************      
!                                                                *      
! SUBROUTINE SUBST solves a linear system of equations.          *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! K       : index of the coordinate function                     *      
! N       : total number of equations and unknowns               *      
! IHF     : 2-dim. array IHF(1:LDIHF,1:N+1); ) auxiliary arrays, *      
! HF      : 2-dim. array   HF(1:LDHF,1:N+3); ) providing storage *      
!                                            ) space             *      
! LDIHF   : leading dimension of IHF as defined in the calling   *      
!           program. LDIHF >= N                                  *      
! LDHF    : leading dimension of HF as defined in the calling    *      
!           program. LDHF >= N                                   *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X1      : N-vector X1(1:N); approximate solution vector        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!                                                                *      
!  source : Brown K.M.: A quadratically convergent Newton-like   *      
!           method based upon Gaussian elimination,              *      
!           Siam J. Numer. Anal., vol 6 (1969), 560 - 569        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Johannes Karfusehr                                 *      
!  editor   : Thomas Eul                                         *      
!  date     : 08.21.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!   declarations                                                        
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER K, N, LDIHF, LDHF 
      DIMENSION IHF (LDIHF, N + 1), HF (LDHF, N + 3), X1 (N) 
!                                                                       
      DO 10 KM = K, 2, - 1 
         KMAX = IHF (KM - 1, N + 1) 
         X1 (KMAX) = 0.0D0 
         DO 20 J = KM, N 
            JSUB = IHF (KM, J) 
            X1 (KMAX) = X1 (KMAX) + HF (KM - 1, JSUB) * X1 (JSUB) 
   20    END DO 
         X1 (KMAX) = X1 (KMAX) + HF (KM - 1, N + 1) 
   10 END DO 
      RETURN 
      END SUBROUTINE SUBST                          

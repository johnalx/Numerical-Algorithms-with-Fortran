![  {Nonlinear Root--Mean--Square Fitting}                              
![  {Nonlinear Root--Mean--Square Fitting}*)                            
      SUBROUTINE SNLFIT (X, Y, W, IWFL, PHI, JNDVT, DVT, MAXIT, PSI,    &
      LDA, M, N, INUM, EPS, A, D, S, C, FV, SQERR, IERR)                
!                                                                       
!*****************************************************************      
!                                                                *      
!  The SUBROUTINE SNLFIT computes a non-linear discrete fitting  *      
!  for M+1 given pairs of values (X(I), Y(I)), I=0, ..., M,      *      
!  possibly using weights, by finding a non-linear approximating *      
!  function defined by N+1 parameters C(K), K=0, ..., M.         *      
!  The model functions have to be provided by the user as a      *      
!  FUNCTION-subroutine. Additionally a SUBROUTINE called DVT may *      
!  be provided by the user (if JNDVT = 0) that determines the    *      
!  partial derivatives with respect to C(K) for K=0, ..., N.     *      
!  Based on an initial approximation for the desired parameters, *      
!  the optimal parameters for the fitting function are           *      
!  determined using the damped Newton-method for non-linear      *      
!  systems. The linear minimization problem that is to be solved *      
!  for each iteration step is solved using Householder           *      
!  transformations.                                              *      
!  Finally, the value of the fitting function is determined for  *      
!  INUM given locations.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X     (N+1)-vector X(0:M) containing the X-values of the nodes*      
!  Y     (N+1)-vector Y(0:M) containing the Y-values at the nodes*      
!  W     (N+1)-vector W(0:M) containing the positive weights     *      
!  IWFL  if IWFL = 0, the nodes will be weighed according to the *      
!        weights in W.                                           *      
!        Otherwise, all weights are set to 1, i.e., there is no  *      
!        need to define values for W                             *      
!  PHI   FUNCTION-subroutine for the model functions defined in  *      
!        the following form:                                     *      
!                                                                *      
!           DOUBLE PRECISION FUNCTION PHI (C, N, X)              *      
!           IMPLICIT DOUBLE PRECISION (A-H,O-Z)                  *      
!           INTEGER N                                            *      
!           DIMENSION C(0:N), X                                  *      
!           --------------------                                 *      
!           PHI = value of the model function at X               *      
!           --------------------                                 *      
!           RETURN                                               *      
!           END                                                  *      
!                                                                *      
!        In the calling program PHI has to be declared as        *      
!        EXTERNAL.                                               *      
!                                                                *      
!  JNDVT if JNDVT = 0, the user has provided a subroutine called *      
!        DVT that determines the partial derivatives.            *      
!        Otherwise the partial derivatives are approximated by   *      
!        central difference quotients.                           *      
!  DVT   user supplied SUBROUTINE that determines the partial    *      
!        derivatives:                                            *      
!                                                                *      
!            SUBROUTINE DVT (X, C, N, F)                         *      
!            IMPLICIT DOUBLE PRECISION (A-H,O-Z)                 *      
!            INTEGER N                                           *      
!            DIMENSION C(0:N), F(0:N)                            *      
!            ----------------------                              *      
!            F(0) = partial derivative w. r. t. C(0) at X        *      
!             .                                                  *      
!             .                                                  *      
!             .                                                  *      
!            F(N) = partial derivative w.r.t C(N) at X           *      
!            ----------------------                              *      
!            RETURN                                              *      
!            END                                                 *      
!                                                                *      
!        In the calling program DVT has to be defined as         *      
!        EXTERNAL.                                               *      
!  MAXIT maximum number of iterations to be performed            *      
!  PSI   INUM-vector PSI(1:INUM) containing the locations where  *      
!        the fitting function is to be evaluated                 *      
!  LDA   leading dimension of the matrix A   (LDA >= M )         *      
!  M     M+1 = number of nodes                                   *      
!  N     N+1 = number of model functions used                    *      
!  INUM  at INUM locations the fitting function is to evaluated  *      
!  EPS   relative error bound for the precision of the optimal   *      
!        parameters                                              *      
!  C     (N+1)-vector C(0:N) containing the initial approxima-   *      
!        tions for the parameters                                *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!                                                                *      
!  A     2-dim array A(0:LDA, 0:N+2), the Jacobi matrix          *      
!  D     (N+1)-vector D(0:N), auxiliary vector for SHOUSE        *      
!  S     (N+1)-vector S(0:N) for the Newton-direction            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  C      (N+1)-vector C(0:N), containing the parameters for the *      
!         fitting function                                       *      
!  FV     INUM-vector FV(1:INUM), containing the values of the   *      
!         fitting function at the INUM given locations           *      
!         PSI(I), I=1(1)INUM                                     *      
!  IERR   error parameter:                                       *      
!         = 0  : everything o.k.                                 *      
!         = 1  : error in the input parameters                   *      
!         = 2  : an error occurred in SUBROUTINE SHOUSE because  *      
!                (numerically) the functional matrix does not    *      
!                have maximal rank                               *      
!         = 3  : after MAXIT iterations the required precision   *      
!                was not reached, i.e., either EPS was chosen    *      
!                too small or the iteration does not converge    *      
!                possibly because the original approximations    *      
!                were too imprecise.                             *      
!  SQERR  least square error                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SNLPRE, SHOUSE, SENORM, MACHPD          *      
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
      DOUBLEPRECISION NEWERR 
      INTEGER IWFL, JNDVT, MAXIT, LDA, M, N, INUM, IERR 
      DIMENSION X (0:M), Y (0:M), W (0:M), PSI (INUM), A (0:LDA, 0:N +  &
      2), D (0:N), S (0:N), C (0:N), FV (INUM)                          
      EXTERNAL DVT, PHI 
      IERR = 0 
!                                                                       
!  testing the input parameters                                         
!                                                                       
      IF (M.LT.N.OR.LDA.LT.M.OR.MAXIT.LE.0.OR.EPS.LE.0.0D0.OR.M.LE.0.OR.&
     &N.LE.0) THEN                                                      
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!  determine the machine constant if the Jacobi matrix                  
!  will be computed via central difference quotients                    
!                                                                       
      IF (JNDVT.NE.0) THEN 
         EPSMA = 1.0D0 
    2    EPSMA = EPSMA * 0.5D0 
         IF (MACHPD (EPSMA + 1.0D0) .EQ.1) GOTO 2 
         EPSMA = EPSMA * 2.0D0 
      ENDIF 
!                                                                       
!  store the weights W in the N+2-nd column of A                        
!                                                                       
      IF (IWFL.EQ.0) THEN 
         DO 5 I = 0, M 
            A (I, N + 2) = DSQRT (W (I) ) 
    5    END DO 
      ELSE 
         DO 6 I = 0, M 
            A (I, N + 2) = 1.0D0 
    6    END DO 
      ENDIF 
!                                                                       
!  form the differences between the Y-values at the nodes               
!  and the values of the fitting function (considering the              
!  weights) and store in the N+1-st column of A                         
!                                                                       
      DO 10 I = 0, M 
         A (I, N + 1) = (Y (I) - PHI (C, N, X (I) ) ) * A (I, N + 2) 
   10 END DO 
!                                                                       
!  determine the square error                                           
!                                                                       
      SQERR = SENORM (A (0, N + 1), M) 
      L = 1 
!                                                                       
!  Newton-iteration:                                                    
!  -----------------                                                    
!  form the functional matrix considering the weights                   
!                                                                       
  100 CALL SNLPRE (X, A (0, N + 2), PHI, DVT, JNDVT, C, LDA, M, N,      &
      EPSMA, D, A)                                                      
!                                                                       
!  determine an improvement S;                                          
!  if in this process an error occurs in SUBROUTINE SHOUSE,             
!  due to the matrix A (numerically) not having maximal rank,           
!  then we set  IERR = 2 and return                                     
!                                                                       
      CALL SHOUSE (A, LDA, M, N, D, S, MARK) 
      IF (MARK.EQ.1) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
!                                                                       
!  damped Newton step                                                   
!                                                                       
      DO 40 K = 0, 10 
         SMOOTH = 1.0D0 / 2.0D0**K 
         DO 20 I = 0, N 
            D (I) = C (I) + S (I) * SMOOTH 
   20    END DO 
         DO 30 I = 0, M 
            A (I, N + 1) = (Y (I) - PHI (D, N, X (I) ) ) * A (I, N + 2) 
   30    END DO 
         NEWERR = SENORM (A (0, N + 1), M) 
         IF (NEWERR.LE.SQERR) GOTO 50 
   40 END DO 
!                                                                       
!  next, the difference between the old and the new                     
!  approximations is stored in S and thereafter the new                 
!  approximation is stored in C                                         
!                                                                       
   50 DO 60 I = 0, N 
         S (I) = C (I) - D (I) 
         C (I) = D (I) 
   60 END DO 
!                                                                       
!  check the stopping criteria                                          
!                                                                       
      IF (SENORM (S, N) .GT.EPS * SENORM (C, N) ) THEN 
         IF (L.LT.MAXIT) THEN 
            SQERR = NEWERR 
            L = L + 1 
            GOTO 100 
         ELSE 
            IERR = 3 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!  if the required precision has been reached, the functional           
!  values of the computed fitting function are computed at              
!  the locations specified in PSI                                       
!                                                                       
      DO 70 I = 1, INUM 
         FV (I) = PHI (C, N, PSI (I) ) 
   70 END DO 
!                                                                       
!  determine least square error                                         
!                                                                       
      SQERR = DSQRT (NEWERR) 
      MAXIT = L 
      RETURN 
      END SUBROUTINE SNLFIT                         

      SUBROUTINE SMNEWT (FX, DFX, N, MAXIT, IERR, KMAX, LUN, EPS,       &
      RNORM2, F, X, DF, LDDF, IWORK, WORK)                              
!                                                                       
!*****************************************************************      
!                                                                *      
!     SMNEWT finds a solution of the nonlinear system of         *      
!     equations                                                  *      
!                   F1(X(1),...,X(N))=0                          *      
!                   F2(X(1),...,X(N))=0                          *      
!                   - - - - - - - - - -                          *      
!                   FN(X(1),...,X(N))=0                          *      
!                                                                *      
!     via the damped Newton method, if it converges for the      *      
!     starting vector.                                           *      
!                                                                *      
!     Three break-off criteria are used:                         *      
!     1.  maximum number of iterations is reached                *      
!     2.  euclidean norm of the difference between the old and   *      
!         new approximate solutions is smaller than or equal to  *      
!         the preset accuracy bound EPS                          *      
!     3.  euclidean norm the function value at the the new       *      
!         approximation is smaller than or equal to EPS          *      
!                                                                *      
!     If desired, output of intermediate results can be          *      
!     generated via input parameter LUN.                         *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! FX      : SUBROUTINE that has to be provided by the user.      *      
!           It defines the system of equations to be solved.     *      
!           In the calling program FX has to be defined as       *      
!           EXTERNAL and must have the form:                     *      
!               SUBROUTINE FX (N,X,F)                            *      
!               DOUBLE PRECISION X(N),F(N)                       *      
!               F(1) = F1 (X(1),...,X(N))                        *      
!               F(2) = F2 (X(1),...,X(N))                        *      
!               - - - - - - - - - - - - - -                      *      
!               F(N) = FN (X(1),...,X(N))                        *      
!               RETURN                                           *      
!               END                                              *      
! DFX     : SUBROUTINE, that has to be provided by the user.     *      
!           It determines the Jacobi matrix of FX.               *      
!           In the calling program DFX has to be defined as      *      
!           EXTERNAL in the form:                                *      
!               SUBROUTINE DFX (N,X,DF,LDDF)                     *      
!               DOUBLE PRECISION DF(LDDF,N), X(N)                *      
!               DF(1,1) = (D F1/D X1) (X(1),...,X(N))            *      
!               .....................................            *      
!               DF(1,N) = (D F1/D XN) (X(1),...,X(N))            *      
!               DF(2,1) = (D F2/D X1) (X(1),...,X(N))            *      
!               .....................................            *      
!               DF(2,N) = (D F2/D XN) (X(1),...,X(N))            *      
!               - - - - - - - - - - - - - - - - - - -            *      
!               DF(LDDF,1)=(D FN/D X1) (X(1),...,X(N))           *      
!               .....................................            *      
!               DF(LDDF,N)=(D FN/D XN) (X(1),...,X(N))           *      
!               RETURN                                           *      
!               END                                              *      
! N       : number of equations and number of unknowns in the    *      
!           given system of equations                            *      
! LUN     : > 0, file number onto which the iteration steps      *      
!                are output                                      *      
!           = 0, no output                                       *      
! MAXIT   : maximum number of iterations to be executed          *      
! KMAX    : damping bound >= 0 ; if KMAX = 0  the standard       *      
!           Newton method is used                                *      
! EPS     : error parameter                                      *      
! X       : N-vector X(1:N); starting vector                     *      
! DF      : 2-dimensional array DF(1:LDDF,1:N); the Jakobi matrix*      
!           (provision of storage space)                         *      
! LDDF    : leading dimension of DF as defined in the calling    *      
!           program. LDDF >= N                                   *      
! IWORK   : N-vector IWORK(1:N); auxiliary vector for the pivot  *      
!           vector for solving the linear system of equations    *      
! WORK    : (4N)-vector WORK(1:4*N); work vector for X, -DELTA X,*      
!           the old X and for the scaling factors in GAUSSP or   *      
!           for the functional values                            *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! MAXIT   : number of iterations executed                        *      
! IERR    : error parameter                                      *      
!           IERR=0 : successful run                              *      
!           IERR=1 : after MAXIT steps the desired accuracy was  *      
!                    not reached                                 *      
!           IERR=2 : error when solving the linear system of     *      
!                    equations (matrix is singular)              *      
!           IERR=3 : incorrect input parameter                   *      
! RNORM2  : accuracy estimate                                    *      
!               RNORM2 = MIN (XNORM2,FNORM2)                     *      
!           (compare description of local variables)             *      
! X       : N-vector X(1:N); approximate solution                *      
! F       : N-vector F(1:N); functional values at the new        *      
!           approximate solution                                 *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! MARK    : error parameter of GAUSSP                            *      
! IX      : starting index for the vector X in WORK              *      
! IDELTA  : starting index for the stepsize DELTA X in WORK      *      
! IXOLD   : starting index for the old approximate solution X in *      
!           WORK                                                 *      
! IFGAUS  : starting index for the scaling factors of GAUSSP in  *      
!           WORK                                                 *      
! IF      : starting index for the functional value vector at X. *      
!           In work this storage space is shared with the one    *      
!           for the scaling factors in GAUSSP.                   *      
! XANRM2  : euclidean norm of F(X)                               *      
! XNNRM2  : euclidean norm of F(X + 1/2**K * DELTA X)            *      
! XNNRMH  : set equal to XNNRM2; however, it is not erased during*      
!           damping. Therefore it can be used later instead of   *      
!           XNNRM2 if there has been no damping; thus XNNRM2 is  *      
!           equal to XNRMH with K=0.                             *      
! XNORM2  : relative accuracy                                    *      
! FNORM2  : euclidean norm of the function at the new approximate*      
!           solution                                             *      
! K       : damping loop counter. The damping factor is 1/2**K   *      
! IT      : Newton iteration loop counter                        *      
! I       : control variable                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GAUSSP, GAUSSS, FENORM                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 05.13.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER MARK, N, MAXIT, IERR, IX, IF, IXOLD, IFGAUS, IT, K, I,    &
      LUN, KMAX, LDDF                                                   
      INTEGER IWORK (N) 
      DOUBLEPRECISION F (N), DF (LDDF, N), X (N), WORK (4 * N) 
!                                                                       
!*****************************************************************      
!*                checking the input parameters                  *      
!*****************************************************************      
!                                                                       
      IF (MAXIT.GT.0.AND.KMAX.GE.0.AND.EPS.GT.0.0D0.AND.N.GT.0.AND.LUN.G&
     &E.0.AND.LDDF.GE.N) THEN                                           
!                                                                       
!*****************************************************************      
!*                        initialization                         *      
!*****************************************************************      
!                                                                       
         IERR = 1 
         IX = 1 
         IDELTA = N + 1 
         IXOLD = 2 * N + 1 
         IFGAUS = 3 * N + 1 
         IF = 3 * N + 1 
         K = 0 
         IT = 0 
!                                                                       
!*****************************************************************      
!*                       Newton iteration                        *      
!*****************************************************************      
!                                                                       
         IF (LUN.GT.0) THEN 
            I = 1 
            WRITE (LUN, 1000) 
            WRITE (LUN, 1100) IT, I, X (1) 
            WRITE (LUN, 1200) (I, X (I), I = 2, N) 
         ENDIF 
!                                                                       
!        calculating the function value and its euclidean               
!        norm at the starting point                                     
!                                                                       
         CALL FX (N, X, F) 
         FNORM2 = FENORM (N, F) 
  100    CONTINUE 
         IT = IT + 1 
!                                                                       
!           calculation of the function value at new approximation.     
!           If damping was used, this amounts only to a relabelling     
!           since the use of norms in the damped algorithm presuppose   
!           knowledge of the new function values                        
!                                                                       
         IF (K.GT.0) THEN 
            DO 10 I = 1, N 
               F (I) = WORK (IF + I - 1) 
   10       END DO 
         ENDIF 
!                                                                       
!           in the new step the euclidean norm of the function value,   
!           FNORM2, at the new approximation becomes the norm of        
!           the function value at the old approximate solution          
!                                                                       
         XANRM2 = FNORM2 
!                                                                       
!           calculation of the Jacobi matrix                            
!                                                                       
         CALL DFX (N, X, DF, LDDF) 
!                                                                       
!*****************************************************************      
!*                  solving  DF * DELTA X = F                    *      
!*****************************************************************      
!                                                                       
!           1. LR factorization                                         
!                                                                       
         CALL GAUSSP (N, DF, LDDF, IWORK, MARK, WORK (IFGAUS) ) 
!                                                                       
!           checking for nonsingularity                                 
!                                                                       
         IF (MARK.NE.0) THEN 
!                                                                       
!              2. solving the linear system of equations by             
!                 using the LR factors from GAUSSP                      
!                                                                       
            CALL GAUSSS (N, DF, LDDF, IWORK, F, WORK (IDELTA) ) 
!                                                                       
!*****************************************************************      
!*  iteration step without damping, saving of the old X, and     *      
!*  determining the euclidean norm of F(X + DELTA X)             *      
!*****************************************************************      
!                                                                       
            DO 20 I = 1, N 
               WORK (IX + I - 1) = X (I) 
               WORK (IXOLD+I - 1) = X (I) 
               X (I) = X (I) - WORK (IDELTA + I - 1) 
   20       END DO 
!                                                                       
!              determination of the norm of F(X + DELTA X)              
!                                                                       
            CALL FX (N, X, F) 
            XNNRM2 = FENORM (N, F) 
            XNNRMH = XNNRM2 
!                                                                       
!*****************************************************************      
!*                           damping                             *      
!*****************************************************************      
!                                                                       
            K = 0 
  200       IF (K.EQ.KMAX.OR.XANRM2.GT.XNNRM2) GOTO 300 
            K = K + 1 
!                                                                       
!                 Newton step with damping                              
!                                                                       
            DO 30 I = 1, N 
               WORK (IDELTA + I - 1) = 0.5D0 * WORK (IDELTA + I - 1) 
               WORK (IX + I - 1) = WORK (IXOLD+I - 1) - WORK (IDELTA +  &
               I - 1)                                                   
   30       END DO 
!                                                                       
!                 determining the euclidean norm of F(X + 1/2**K * DELTA
!                                                                       
            CALL FX (N, WORK (IX), WORK (IF) ) 
            XNNRM2 = FENORM (N, WORK (IF) ) 
            GOTO 200 
  300       CONTINUE 
!                                                                       
!              if XANRM2 > XNNRM2, calculations are continued           
!              with damping; K=0 indicates that there was               
!              no damping                                               
!                                                                       
            IF (XANRM2.GT.XNNRM2.AND.K.GT.0) THEN 
               DO 40 I = 1, N 
                  X (I) = WORK (IX + I - 1) 
   40          END DO 
            ELSE 
               K = 0 
            ENDIF 
!                                                                       
!*****************************************************************      
!*          test for accuracy and possible stop                  *      
!*****************************************************************      
!                                                                       
!              1. checking the second break-off criterion,              
!                 i.e., determining XNORM2                              
!                                                                       
            DO 50 I = 1, N 
               WORK (IX + I - 1) = X (I) - WORK (IXOLD+I - 1) 
   50       END DO 
            XNORM2 = FENORM (N, X) 
            IF (XNORM2.GT.0) THEN 
               XNORM2 = FENORM (N, WORK (IX) ) / XNORM2 
            ELSE 
               XNORM2 = FENORM (N, WORK (IX) ) 
            ENDIF 
!                                                                       
!              2. checking the third break-off criterion,               
!                 i.e., determining FNORM2                              
!                                                                       
            IF (K.GT.0) THEN 
               FNORM2 = XNNRM2 
            ELSE 
               FNORM2 = XNNRMH 
            ENDIF 
!                                                                       
!              the smallest of the two error estimates is chosen        
!                                                                       
            IF (XNORM2.GT.FNORM2) THEN 
               RNORM2 = FNORM2 
            ELSE 
               RNORM2 = XNORM2 
            ENDIF 
            IF (LUN.GT.0) THEN 
               I = 1 
               WRITE (LUN, 1100) IT, I, X (1), RNORM2, K 
               WRITE (LUN, 1200) (I, X (I), I = 2, N) 
            ENDIF 
            IF (RNORM2.LE.EPS) THEN 
               MAXIT = IT 
               IERR = 0 
            ENDIF 
         ELSE 
            IERR = 2 
            MAXIT = IT 
         ENDIF 
         IF (IT.LT.MAXIT.AND.IERR.EQ.1) GOTO 100 
      ELSE 
         IERR = 3 
         MAXIT = 0 
      ENDIF 
!                                                                       
!*****************************************************************      
!*                           formats                             *      
!*****************************************************************      
!                                                                       
 1000 FORMAT (1X,'ITERATION STEP',10X,'APROXIMATION',14X,               &
     &           'ACCURACY ESTIMATE   K')                               
 1100 FORMAT (1X,I6,2X,I6,8X,D22.15,6X,D22.15,2X,I3) 
 1200 FORMAT (1X, 6X,2X,I6,8X,D22.15) 
!                                                                       
      RETURN 
      END SUBROUTINE SMNEWT                         

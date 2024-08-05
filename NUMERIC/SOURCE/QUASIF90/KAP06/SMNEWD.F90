![KA{P 6}{Systems of Nonlinear Equations}                               
![       {Systems of Nonlinear Equations}*)                             
![             {Damped Newton Method for Systems}*)                     
      SUBROUTINE SMNEWD (FX, N, MAXIT, IERR, KMAX, LUN, IUPD, EPS,      &
      RNORM2, F, X, DF, LDDF, IWORK, WORK)                              
!                                                                       
!*****************************************************************      
!                                                                *      
!     SMNEWD finds a solution of a nonlinear system of equations *      
!                                                                *      
!                   F1(X(1),...,X(N))=0                          *      
!                   F2(X(1),...,X(N))=0                          *      
!                   - - - - - - - - - -                          *      
!                   FN(X(1),...,X(N))=0                          *      
!                                                                *      
!     by the damped Newton method, if this method converges for  *      
!     the starting vector. Here the Jacobi matrix is replaced    *      
!     by the forward difference quotients.                       *      
!     Hence the user does not need to supply the partial         *      
!     derivatives. The parameter IUPD can be used to predetermine*      
!     after how many iterations the Jakobi matrix is to be re-   *      
!     computed and LR decomposed.                                *      
!     IUPD = 1 specifies the damped Newton method.               *      
!     IUPD > 1 specifies a damped version of the simplified      *      
!              Newton method.                                    *      
!     In general, more iteration steps are required if IUPD > 1. *      
!     However, since the function then does not have to be       *      
!     evaluated and the Jacobi matrix not LR factored as often,  *      
!     this may save computational time.                          *      
!                                                                *      
!     Three break-off criteria are used:                         *      
!     1.  maximum number of iterations has been reached          *      
!     2.  euclidean norm of the difference between the old and   *      
!         the new approximate solutions is smaller or equal to   *      
!         the preset (relative) accuracy bound EPS               *      
!     3.  euclidean norm of the function value at the new        *      
!         approximate solution is smaller or equal to EPS        *      
!                                                                *      
!     If desired, there will be output of intermediate results   *      
!     (via input parameter LUN).                                 *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! FX      : SUBROUTINE, that has to be specified by the user. It *      
!           describes the system of equations to be solved.      *      
!           In the calling program FX has to be defined as       *      
!           EXTERNAL. It has to be in the form:                  *      
!               SUBROUTINE FX (N,X,F)                            *      
!               DOUBLE PRECISION X(N),F(N)                       *      
!               F(1) = F1 (X(1),...,X(N))                        *      
!               F(2) = F2 (X(1),...,X(N))                        *      
!               - - - - - - - - - - - - -                        *      
!               F(N) = FN (X(1),...,X(N))                        *      
!               RETURN                                           *      
!               END                                              *      
! N       : number of equations and number of unknowns           *      
! LUN     : > 0, file number onto which the iterates are stored  *      
!           = 0, no output                                       *      
! MAXIT   : maximum number of iterations to be executed          *      
! KMAX    : a bound for the damping factor >= 0;                 *      
!           KMAX = 0  ==> standard Newton method is used,        *      
!           usually a value between 4 and 6 is chosen for KMAX   *      
! IUPD    : after each IUPD steps the Jacobi matrix is           *      
!           reconfigured and LR decomposed anew. In general, the *      
!           method will not converge if IUPD is chosen too large.*      
!           IUPD between 1 and 4 usually will give meaningful    *      
!           results.                                             *      
! EPS     : error parameter                                      *      
! X       : N-vector X(1:N); starting vector                     *      
! DF      : 2-dim. array DF(1:LDDF,1:N); Jacobi matrix           *      
!           (providing storage space)                            *      
! LDDF    : leading dimension of DF as defined in the calling    *      
!           program. LDDF >= N                                   *      
! IWORK   : N-vector IWORK(1:N); auxiliary vector, the pivot     *      
!           vector for solving the linear system of equations    *      
! WORK    : (4N)-vector WORK(1:4*N); auxiliary vector for X,     *      
!           -DELTA X, the old iterate X and the scaling factors  *      
!           in GAUSSP or the functional values at X              *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! MAXIT   : number of iterations executed                        *      
! IERR    : error parameter                                      *      
!           IERR=0 : no error                                    *      
!           IERR=1 : error bound was not reached after MAXIT     *      
!                    steps                                       *      
!           IERR=2 : error when solving the linear system of     *      
!                    equations (matrix singular)                 *      
!           IERR=3 : incorrect input parameter                   *      
! RNORM2  : accuracy estimate                                    *      
!               RNORM2 = MIN (XNORM2,FNORM2)                     *      
!           (compare with the description of the local variables)*      
! X       : N-vector X(1:N); approximate solution                *      
! F       : N-vector F(1:N); function value at the approximate   *      
!           solution                                             *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! MARK    : error parameter from GAUSSP                          *      
! IX      : starting index for X in the auxiliary vector WORK    *      
! IDELTA  : starting index for the stepsize DELTA X in WORK      *      
! IXOLD   : starting index for the old approximate solution X    *      
!           in WORK                                              *      
! IFGAUS  : starting index for the scaling factors from GAUSSP   *      
!           in WORK                                              *      
! IF      : starting index for the vector of functional values   *      
!           at X. In WORK this storage space is identical and    *      
!           shared with the one for the scaling factors in       *      
!           GAUSSP.                                              *      
! XANRM2  : euclidean norm of F(X)                               *      
! XNNRM2  : euclidean norm of F(X + 1/2**K * DELTA X)            *      
! XNNRMH  : set equal to XNNRM2; however, this constant it is not*      
!           altered during damping and thus it can be used later *      
!           instead of XNNRM2 if damping will not be used. Thus  *      
!           it equals XNNRM2 if K=0.                             *      
! XNORM2  : relative accuracy                                    *      
! FNORM2  : euclidean norm of the function value at the new      *      
!           approximation                                        *      
! K       : damping loop counter. The damping factor is 1/2**K   *      
! IT      : Newton-iteration loop counter                        *      
! I       : control variable                                     *      
! IUP     : counter for estimating the Jacobi matrix             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GAUSSP, GAUSSS, FENORM, FDIFQU, MACHPD  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 07.02.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER IWORK (N) 
      DOUBLEPRECISION F (N), DF (LDDF, N), X (N), WORK (4 * N) 
      EXTERNAL FX 
!                                                                       
!*****************************************************************      
!*               checking the input parameters                   *      
!*****************************************************************      
!                                                                       
      IF (MAXIT.GT.0.AND.KMAX.GE.0.AND.EPS.GT.0.0D0.AND.N.GT.0.AND.LUN.G&
     &E.0.AND.IUPD.GT.0.AND.LDDF.GE.N) THEN                             
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
         IUP = IUPD-1 
         K = 0 
         IT = 0 
!                                                                       
!*****************************************************************      
!*                      Newton iteration                         *      
!*****************************************************************      
!                                                                       
         IF (LUN.GT.0) THEN 
            I = 1 
            WRITE (LUN, 1000) 
            WRITE (LUN, 1100) IT, I, X (1) 
            WRITE (LUN, 1200) (I, X (I), I = 2, N) 
         ENDIF 
!                                                                       
!        determining the functional value and its euclidean norm        
!        at the starting point                                          
!                                                                       
         CALL FX (N, X, F) 
         FNORM2 = FENORM (N, F) 
  100    CONTINUE 
         IT = IT + 1 
!                                                                       
!           determining the functional value at the new                 
!           approximate solution. However, this only involves           
!           relabelling if we have used damping, since the norms        
!           necessary for the damped algorithm automatically            
!           presuppose knowledge of the functional value                
!                                                                       
         IF (K.GT.0) THEN 
            DO 10 I = 1, N 
               F (I) = WORK (IF + I - 1) 
   10       END DO 
         ENDIF 
!                                                                       
!           the euclidean norm of the function value at the new         
!           approximation, FNORM2, becomes the euclidean norm of        
!           the functional value for the old approximate solution       
!                                                                       
         XANRM2 = FNORM2 
!                                                                       
!*****************************************************************      
!*                  solving of DF * DELTA X = F                  *      
!*****************************************************************      
!                                                                       
!           if necessary estimation and multiplication                  
!           of the Jacobi matrix                                        
!                                                                       
         IUP = IUP + 1 
         IF (IUP.EQ.IUPD) THEN 
            IUP = 0 
            CALL FDIFQU (FX, N, X, F, WORK (IF), DF, LDDF) 
!                                                                       
!              1. LR factorization                                      
!                                                                       
            CALL GAUSSP (N, DF, LDDF, IWORK, MARK, WORK (IFGAUS) ) 
         ENDIF 
!                                                                       
!           checking for singularity                                    
!                                                                       
         IF (MARK.NE.0) THEN 
!                                                                       
!              2. solving the linear system of equations with           
!                 the LR factors from GAUSSP                            
!                                                                       
            CALL GAUSSS (N, DF, LDDF, IWORK, F, WORK (IDELTA) ) 
!                                                                       
!*****************************************************************      
!*  iteration step without damping, saving of the old X, and     *      
!*  determining the norm of F(X + DELTA X)                       *      
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
!*                          damping                              *      
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
!                 determining the norm of F(X + 1/2**K * DELTA X)       
!                                                                       
            CALL FX (N, WORK (IX), WORK (IF) ) 
            XNNRM2 = FENORM (N, WORK (IF) ) 
            GOTO 200 
  300       CONTINUE 
!                                                                       
!              if XANRM2 > XNNRM2, computations are continued           
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
!*         test for accuracy and if warranted : stop             *      
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
!              we make use of the smallest of the two error estimates   
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
!*                    output formats                             *      
!*****************************************************************      
!                                                                       
 1000 FORMAT (1X,'ITERATION STEP',10X,'APROXIMATION',14X,               &
     &           'ERROR ESTIMATE  K')                                   
 1100 FORMAT (1X,I6,2X,I6,8X,D22.15,6X,D22.15,2X,I3) 
 1200 FORMAT (1X, 6X,2X,I6,8X,D22.15) 
!                                                                       
      RETURN 
      END SUBROUTINE SMNEWD                         
!                                                                       
!                                                                       
      SUBROUTINE FDIFQU (FX, N, X, F, FXJPH, DF, LDDF) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Approximates the  Jacobi matrix for the function               *      
!                                                                *      
!                  F1 (X(1),...,X(N))                            *      
!                  F2 (X(1),...,X(N))                            *      
!                  - - - - - - - - - -                           *      
!                  FN (X(1),...,X(N))                            *      
!                                                                *      
! by forward difference quotients.                               *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! FX      : SUBROUTINE, that has to be provided by the user.     *      
!           It finds the functional value at X. The program      *      
!           computes an approximate 1st derivative of FX at X.   *      
!           In the calling program FX has to be defined as       *      
!           EXTERNAL of the form:                                *      
!               SUBROUTINE FX (N,X,F)                            *      
!               DOUBLE PRECISION X(N),F(N)                       *      
!               F(1) = F1 (X(1),...,X(N))                        *      
!               F(2) = F2 (X(1),...,X(N))                        *      
!               - - - - - - - - - - - - -                        *      
!               F(N) = FN (X(1),...,X(N))                        *      
!               RETURN                                           *      
!               END                                              *      
! N       : number of component functions and variables of FX    *      
! X       : N-vector X(1:N); location where the 1st derivative   *      
!           is to be estimated                                   *      
! F       : N-vector F(1:N); functional value of FX at X         *      
! FXJPH   : N-vector FXJPH(1:N); provides storage space for the  *      
!           function values                                             
!               F1(X(1),...,X(J-1),X(J)+H,X(J+1),...,X(N))       *      
!               F2(X(1),...,X(J-1),X(J)+H,X(J+1),...,X(N))       *      
!               - - - - - - - - - - - - - - - - - - - - - -      *      
!               FN(X(1),...,X(J-1),X(J)+H,X(J+1),...,X(N))       *      
! LDDF    : leading dimension of DF as defined in the calling    *      
!           program. LDDF >= N                                   *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! DF      : 2-dimensional array DF(1:LDDF,1:N); the computed     *      
!           Jacobi matrix at X                                   *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! H       : stepsize for the forward difference quotients, H is  *      
!           taken as the square root of the machine constant     *      
! HBEST   : logical variable for determining H. Initially HBEST  *      
!           is set to .TRUE., so that on the initial call of     *      
!           FDIFQU the constant H is determined. Then HBEST is   *      
!           set to .FALSE., and H and HBEST are stored unchanged *      
!           by further calls of FDIFQU                           *      
! HAEPSM  : auxiliary variable for determining the machine       *      
!           constant                                             *      
! EPSM    : machine constant                                     *      
! XJ      : auxiliary variable used for forming the difference   *      
!           quotient with respect to X(J)                        *      
! I       : control variable                                     *      
! J       : control variable                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 02.07.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!     declarations                                                      
!                                                                       
      INTEGER N, LDDF 
      DOUBLEPRECISION X (N), F (N), FXJPH (N), DF (LDDF, N) 
!                                                                       
!     local variables                                                   
!                                                                       
      INTEGER I, J 
      LOGICAL HBEST 
!                                                                       
      SAVE H, HBEST 
!                                                                       
      DATA HBEST / .TRUE. / 
!                                                                       
      IF (HBEST) THEN 
!                                                                       
!*****************************************************************      
!*          determining the square root of the machine constant  *      
!*****************************************************************      
!                                                                       
!        EPSM represents the smallest positive number for which         
!        (1.0+EPSM) .GT. 1.0 . EPSM is determined as a power of 1./2.   
!                                                                       
         HAEPSM = 1.0D0 
   10    CONTINUE 
         HAEPSM = 0.5D0 * HAEPSM 
         IF (MACHPD (1.0D0 + HAEPSM) .EQ.1) GOTO 10 
         EPSM = 2.0D0 * HAEPSM 
         H = DSQRT (EPSM) 
         HBEST = .FALSE. 
      ENDIF 
!                                                                       
      DO 30 J = 1, N 
         XJ = X (J) 
         X (J) = XJ + H 
         CALL FX (N, X, FXJPH) 
         X (J) = XJ 
         DO 20 I = 1, N 
            DF (I, J) = (FXJPH (I) - F (I) ) / H 
   20    END DO 
   30 END DO 
      RETURN 
      END SUBROUTINE FDIFQU                         

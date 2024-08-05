      SUBROUTINE NLESYS (FX, M, N, JACMA, DFX, DF, LDDF, MAXIT, EPS,    &
      KMAX, LUN, X, F, RNORM2, IERR, D, S, WORK)                        
!                                                                       
!*****************************************************************      
!                                                                *      
!  SUBROUTINE NLESYS solves a nonlinear system of equations      *      
!  consisting of M+1 equations in N+1 unknowns by applying the   *      
!  damped Newton method provided that it converges for the       *      
!  starting vector.                                              *      
!  The Jacobi matrix must be either determined by a subroutine   *      
!  provided by the user or approximated by central difference    *      
!  quotients.                                                    *      
!  Three break-off criteria are used:                            *      
!    1. maximum number of iterations has been reached            *      
!    2. the relative difference between the euclidean norms of   *      
!       the old and new approximations is smaller or equal to    *      
!       EPS                                                      *      
!    3. the euclidean norm of the function value for the newest  *      
!       approximation is smaller than or equal to EPS            *      
!                                                                *      
!  this program was developed by using the SUBROUTINES SMNEWD    *      
!  and SMNEWT by Thomas Eul and by SUBROUTINE SNLFIT by          *      
!  Ilona Westermann.                                             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  ================+                                             *      
!  FX     : Function SUBROUTINE, that has to be provided by the  *      
!           user. It defines the system of equations that are to *      
!           be solved. In the calling program FX has to be       *      
!           defined as EXTERNAL. It must have the form:          *      
!               SUBROUTINE FX (X, N, F, M)                       *      
!               INTEGER M, N                                     *      
!               DOUBLE PRECISION X(0:N), F(0:M)                  *      
!               ------------------                               *      
!               F(0) = F0 (X(0), ... , X(N)                      *      
!                .                                               *      
!                .                                               *      
!               F(M) = FM (X(0), ... , X(N)                      *      
!               ------------------                               *      
!               RETURN                                           *      
!               END                                              *      
!  M      : M+1 = number of equations                            *      
!  N      : N+1 = number of unknowns                             *      
!  JACMA  : logical parameter:                                   *      
!           If JACMA = .TRUE., the user has provided his own     *      
!           subroutine for determining the Jacobi matrix.        *      
!           Otherwise the Jacobi matrix is computed intrinsically*      
!           by central difference quotients.                     *      
!  DFX    : A SUBROUTINE that has to be provided by the user if  *      
!           JACMA = .TRUE.  DFX determines the Jacobi matrix.    *      
!           In the calling program DFX has to be defined as      *      
!           EXTERNAL and should have the following form:         *      
!               SUBROUTINE DFX (X, M, N, DF, LDDF)               *      
!               INTEGER M, N, LDDF                               *      
!               DOUBLE PRECISION X(0:N), DF(0:LDDF,0:N)          *      
!               ----------------------                           *      
!               DF(I,K) = partial derivative of function I       *      
!                         by X(K) at the position X              *      
!                         K=0(1)N, I=0(1)M                       *      
!               ----------------------                           *      
!               RETURN                                           *      
!               END                                              *      
!  DF     : 2-dim. array DF(0:LDDF,0:N+1); used for storage      *      
!           space for the Jacobi matrix                          *      
!  LDDF   : leading dimension of DF as defined in the calling    *      
!           program                                              *      
!  X      : (N+1)-vector X(0:N) containing the starting vector   *      
!  MAXIT  : maximum number of iterations to be executed          *      
!  EPS    : relative error bound                                 *      
!  KMAX   : damping bound, if KMAX=0 => standard Newton method   *      
!  LUN    : > 0, file number onto which the iteration steps are  *      
!                stored                                          *      
!           = 0, no output                                       *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  D      : (N+1)-vector D(0:N); auxiliary vector for SUBROUTINE *      
!                                SHOUSE                          *      
!  S      : (N+1)-vector S(0:N); for the Newton direction        *      
!  WORK   : (N+1)-vector WORK(0:M); auxiliary vector for         *      
!                                   SUBROUTINE JACOBI            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X      : (N+1)-vector X(0:N) containing the approximate       *      
!           solution                                             *      
!  F      : (N+1)-vector F(0:N) containing the function values   *      
!           at the approximate solution                          *      
!  IERR   : error parameter:                                     *      
!           = 0  : everything o.k.                               *      
!           = 1  : error in the input parameters                 *      
!           = 2  : error while solving the linear minimization   *      
!                  problem (matrix singular)                     *      
!           = 3  : error bound not attained after MAXIT          *      
!                  iterations                                    *      
!  RNORM2 : an estimate for the achieved accuracy,               *      
!           RNORM2 = MIN (XNORM2, FNORM2)                        *      
!           (compare with the description of local variables)    *      
!                                                                *      
!                                                                *      
!  LOCAL VARIABLES:                                              *      
!  ================                                              *      
!  FNORM2 : euclidean norm of F(X)                               *      
!  FNORMN : euclidean norm of F(X + 1/2**K * DELTA X)            *      
!  XNORM2 : relative accuracy                                    *      
!  EPSMA  : machine constant                                     *      
!  SMOOTH : damping factor 1/2**K                                *      
!  MARK   : error parameter of SHOUSE                            *      
!  K      : counter for the damping loop. The damping factor is  *      
!           1/2**K                                               *      
!  IT     : Newton iteration loop counter                        *      
!  I      : control variable                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: JACOBI, SHOUSE, SENORM, MACHPD          *      
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
      INTEGER M, N, LDDF, MAXIT, KMAX, LUN, IERR 
      LOGICAL JACMA 
      DIMENSION DF (0:LDDF, 0:N + 1), X (0:N), F (0:M), D (0:N),        &
      S (0:N), WORK (0:M)                                               
      EXTERNAL FX 
!                                                                       
!  test of the input parameters                                         
!                                                                       
      IF (M.GE.N.AND.LDDF.GE.M.AND.MAXIT.GT.0.AND.KMAX.GE.0.AND.EPS.GT.0&
     &.0D0.AND.N.GE.0.AND.M.GE.0.AND.LUN.GE.0) THEN                     
         IERR = 3 
         IT = 0 
!                                                                       
!  determine machine constant in case the Jacobi matrix is calculated in
!  this program by central difference quotients                         
!                                                                       
         IF (.NOT.JACMA) THEN 
            EPSMA = 1.0D0 
   10       EPSMA = EPSMA * 0.5D0 
            IF (MACHPD (EPSMA + 1.0D0) .EQ.1) GOTO 10 
            EPSMA = EPSMA * 2.0D0 
         ENDIF 
         IF (LUN.GT.0) THEN 
            I = 0 
            WRITE (LUN, 1000) 
            WRITE (LUN, 1100) IT, I, X (0) 
            WRITE (LUN, 1200) (I, X (I), I = 1, N) 
            WRITE (LUN, 1300) 
         ENDIF 
!                                                                       
!  determine functional value and its norm at the starting point        
!                                                                       
         CALL FX (X, N, F, M) 
         FNORM2 = DSQRT (SENORM (F, M) ) 
!                                                                       
!  Newton- teration                                                     
!                                                                       
  100    IT = IT + 1 
!                                                                       
!  determining the Jacobi matrix                                        
!                                                                       
         IF (JACMA) THEN 
            CALL DFX (X, M, N, DF, LDDF) 
         ELSE 
            CALL JACOBI (FX, X, M, N, DF, LDDF, EPSMA, WORK) 
         ENDIF 
!                                                                       
!  minimize (DF * DELTA X + F)      (DELTA X  here is called S)         
!                                                                       
         DO 15 I = 0, M 
            DF (I, N + 1) = - F (I) 
   15    END DO 
         CALL SHOUSE (DF, LDDF, M, N, D, S, MARK) 
!                                                                       
!  check for nonsingularity                                             
!                                                                       
         IF (MARK.NE.1) THEN 
!                                                                       
!  damped Newton step                                                   
!                                                                       
            DO 40 K = 0, KMAX 
               SMOOTH = 1.0D0 / 2.0D0**K 
               DO 20 I = 0, N 
                  D (I) = X (I) + S (I) * SMOOTH 
   20          END DO 
               CALL FX (D, N, F, M) 
               FNORMN = DSQRT (SENORM (F, M) ) 
!                                                                       
!  if FNORMN < FNORM2 continue calculations with damping                
!                                                                       
               IF (FNORMN.LT.FNORM2) THEN 
                  DO 30 I = 0, N 
                     S (I) = S (I) * SMOOTH 
                     X (I) = D (I) 
   30             END DO 
                  GOTO 60 
               ENDIF 
   40       END DO 
!                                                                       
!  if the loop was completed, i.e., damping gave no improvement,        
!  continue calculations without damping                                
!                                                                       
            DO 50 I = 0, N 
               X (I) = X (I) + S (I) 
   50       END DO 
            CALL FX (X, N, F, M) 
            FNORMN = DSQRT (SENORM (F, M) ) 
            K = 0 
   60       CONTINUE 
!                                                                       
!  check accuracy and if necessary stop:                                
!                                                                       
!  1. use relative error estimate for X                                 
!                                                                       
            XNORM2 = DSQRT (SENORM (X, N) ) 
            IF (XNORM2.EQ.0.0D0) THEN 
               XNORM2 = DSQRT (SENORM (S, N) ) 
            ELSE 
               XNORM2 = DSQRT (SENORM (S, N) ) / XNORM2 
            ENDIF 
!                                                                       
!  2. find the norm of the functional value for X                       
!                                                                       
            FNORM2 = FNORMN 
!                                                                       
!     select the smaller one of both norms                              
!                                                                       
            RNORM2 = DMIN1 (FNORM2, XNORM2) 
            IF (LUN.GT.0) THEN 
               I = 0 
               WRITE (LUN, 1100) IT, I, X (0), RNORM2, K 
               WRITE (LUN, 1200) (I, X (I), I = 1, N) 
               WRITE (LUN, 1300) 
            ENDIF 
            IF (RNORM2.LE.EPS) THEN 
               IERR = 0 
               MAXIT = IT 
            ENDIF 
         ELSE 
            IERR = 2 
            MAXIT = IT 
         ENDIF 
         IF (IT.LT.MAXIT.AND.IERR.EQ.3) GOTO 100 
      ELSE 
         IERR = 1 
         MAXIT = 0 
      ENDIF 
      RETURN 
 1000 FORMAT (1X,'ITERATION NUMBER',10X,'APROXIMATION',14X,             &
     &           'ACCURACY ESTIMATE   K')                               
 1100 FORMAT (1X,I6,2X,I6,8X,D22.15,6X,D22.15,2X,I3) 
 1200 FORMAT (1X, 6X,2X,I6,8X,D22.15) 
 1300 FORMAT (1X) 
      END SUBROUTINE NLESYS                         

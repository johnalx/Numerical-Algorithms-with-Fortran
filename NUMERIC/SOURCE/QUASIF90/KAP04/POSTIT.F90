      SUBROUTINE POSTIT (N, A0, A, LDA, IPIVOT, Y, X, EPS, MAXIT, NUMIT,&
      IERR, Z, R, RS)                                                   
!                                                                       
!*****************************************************************      
!                                                                *      
!  POSTIT performs iterative refinement after GAUSSP and GAUSSS  *      
!  have been executed. The iteration is stopped if a set maximum *      
!  MAXIT of iteration steps has been performed or if the relative*      
!  improvement satisfies  MAX(ABS(Z(I)))/MAX(ABS(X(I))) < EPS    *      
!  for I=1, ..., N.                                              *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : order of the matrices A and A0.                    *      
!  A0       : 2-dimensional array A0(1:LDA,1:N); the matrix      *      
!             A(ORG).                                            *      
!  A        : 2-dimensional array A(1:LDA,1:N) containing the    *      
!             factors  L  and  R  with  P * A(ORG) = L * R.      *      
!             P = permutation matrix. A is an output of          *      
!             SUBROUTINE GAUSSP.                                 *      
!  LDA      : leading dimension of A and A0 as defined in the    *      
!             calling program.                                   *      
!  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row permuta-*      
!             tions in  P * A(ORG)  compared with A(ORG). It is  *      
!             an output vector of SUBROUTINE GAUSSP.             *      
!  Y        : N-vector Y(1:N); the right side of the system of   *      
!             equations.                                         *      
!  X        : N-vector X(1:N); the solution of the system of     *      
!             equations; X serves as starting vector for the     *      
!             iterative refinement.                              *      
!  EPS      : error bound for the relative improvement;          *      
!             if EPS < 4 * machine constant, the program inter-  *      
!             ally sets EPS = 4 * machine constant.              *      
!  MAXIT    : maximum number of iterations.                      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X        : N-vector X(1:N); the solution of the linear system.*      
!  EPS      : error bound actually used.                         *      
!  NUMIT    : number of iteration steps executed.                *      
!  IERR     : error parameter.                                   *      
!             = 0 : program stopped since the relative           *      
!                   improvement is < EPS.                        *      
!             = 1 : the set accuracy was not achieved after      *      
!                   MAXIT iterations.                            *      
!             = 2 : same as IERR=1, except that the components of*      
!                   the correction vector have begone to increase*      
!                   indicating divergence and ill-conditioning of*      
!                   A0.                                          *      
!                                                                *      
!                                                                *      
!  AUXILIARY VECTORS:                                            *      
!  ==================                                            *      
!  Z        : N-vector Z(1:N).                                   *      
!  R        : N-vector R(1:N) in DOUBLE PRECISION.               *      
!  RS       : N-vector RS(1:N).                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GAUSSS, MACHPD                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors  : Gisela Engeln-Muellges, Guido Dubois               *      
!  date     : 04.25.88                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A0 (1:LDA, 1:N), A (1:LDA, 1:N), Y (1:N), X (1:N),&
      Z (1:N), RS (1:N)                                                 
      DOUBLEPRECISION R (1:N) 
      INTEGER IPIVOT (1:N) 
!                                                                       
!  local storage of the minimally acceptable error bound EPSMIN         
!  in case that the subroutine is called more than once.                
!                                                                       
      SAVE EPSMIN, IFLAG 
      DATA IFLAG / 0 / 
!                                                                       
!  calculating the machine constant and initializing the minimally      
!  acceptable error bound EPSMIN for the relative improvement.          
!                                                                       
      IF (IFLAG.EQ.0) THEN 
         IFLAG = 1 
         FMACHP = 1.0D0 
   10    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
         EPSMIN = 8.0D0 * FMACHP 
      ENDIF 
      IF (EPS.LT.EPSMIN) EPS = EPSMIN 
      NUMIT = 0 
      ZMA = 1.0D20 
      IERR = 0 
   50 NUMIT = NUMIT + 1 
!                                                                       
!  calculation of the residual vector.                                  
!                                                                       
      DO 20 I = 1, N 
         R (I) = Y (I) 
         DO 30 K = 1, N 
            R (I) = R (I) - A0 (I, K) * X (K) 
   30    END DO 
         RS (I) = R (I) 
   20 END DO 
!                                                                       
!  calculating the correction vector and improving                      
!  the approximate solution.                                            
!                                                                       
      CALL GAUSSS (N, A, LDA, IPIVOT, RS, Z) 
      XM = 0.0D0 
      ZM = 0.0D0 
      DO 40 I = 1, N 
         XM = DMAX1 (XM, DABS (X (I) ) ) 
         ZM = DMAX1 (ZM, DABS (Z (I) ) ) 
         X (I) = X (I) + Z (I) 
   40 END DO 
      IF (ZM / XM.LT.EPS) RETURN 
      IF (NUMIT.GE.MAXIT) THEN 
         IERR = 1 
         IF (ZM.GT.ZMA) IERR = 2 
         RETURN 
      ENDIF 
      ZMA = ZM 
      GOTO 50 
      END SUBROUTINE POSTIT                         

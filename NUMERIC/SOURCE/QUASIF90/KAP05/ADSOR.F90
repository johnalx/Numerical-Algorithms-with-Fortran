![KA{P 5}{Iterative Methods for Linear Systems}                         
![       {Iterative Methods for Linear Systems}*)                       
![  {The Gau"s--Seidel Iteration}{The Gau"s--Seidel Iteration}*)        
      SUBROUTINE ADSOR (A, N, IA, B, X, KADAPT, EPS, KMAX, IMETH,       &
      ISWITC, OMEGA, WORK, RES, ITNUMB, IERR)                           
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program solves an inhomogeneous linear system AX = B of  *      
!  equations with a nonsingular system matrix A. The method of   *      
!  Jacobi is used jointly with relaxation, where the relaxation  *      
!  parameter OMEGA is adjusted during the iteration (adaptive    *      
!  SOR method).                                                  *      
!  For a suitable choice of parameters (refer to the remark      *      
!  below), this program can perform the Gau·-Seidel method or    *      
!  a non-adaptive SOR method.                                    *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : 2-dimensional array A(1:IA,1:N), containing the      *      
!           system matrix for the linear equations               *      
!  N      : size of the linear system                            *      
!  IA     : leading dimension of A, as specified in the calling  *      
!           program                                              *      
!  B      : N-vector B(1:N), the right hand side of the system   *      
!  X      : N-vector X(1:N) containing the starting value for    *      
!           iteration                                            *      
!  KADAPT : Number of iterations, after which the relaxation     *      
!           parameter is to be redefined                         *      
!  EPS    : desired accuracy; the iteration is stopped when the  *      
!           maximum norm of the relative error does not exceed   *      
!           EPS                                                  *      
!  KMAX   : Maximal number of iterations allowed                 *      
!  IMETH  : parameter that determines the method used:           *      
!           = 0, adaptive SOR method                             *      
!           = 1, SOR method for a given relaxation parameter     *      
!           = 2, Gau·-Seidel method                              *      
!  ISWITC : parameter that determines the convergence criterion  *      
!           to be used:                                          *      
!           = 0, none                                            *      
!           = 1, row sum criterion                               *      
!           = 2, column sum criterion                            *      
!           = 3, criterion of Schmidt and v. Mises               *      
!  OMEGA  : in case IMETH=1, the optimal relaxation parameter    *      
!           must be part of the input; otherwise only the name   *      
!           must be declared in the callimng program.            *      
!                                                                *      
!                                                                *      
!  REMARKS:                                                      *      
!  ========                                                      *      
!  For the adaptive SOR method (IMETH=0) we recommend to set     *      
!  KADAPT=4 or KADAPT=5.                                         *      
!  If the optimal relaxationcoefficient Wopt is known for A, then*      
!  one should set IMETH=1 and OMEGA = Wopt, i.e., the SOR method *      
!  with given optimal relaxation coefficient should be used.     *      
!  If IMETH=2, then the program performs the Gau·-Seidel method. *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  WORK   : 2-dim. array  WORK(1:N,1:3)                          *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A      : 2-dim. array A(1:IA,1:N), the input matrix A is over-*      
!           written by: A(I,J)=A(I,J)/A(I,I) for I,J=1, ..., N   *      
!  B      : N-vector B(1:N), the right hand side is replaced by  *      
!           B(I)=B(I)/A(I,I); I=1,N                              *      
!  OMEGA  : - if IMETH = 0, the program returns the adaptively   *      
!             computed relaxations parameter.                    *      
!           - if IMETH = 1, the optimal relaxation parameter     *      
!             is returned as put in externally.                  *      
!           - if IMETH = 2, then on output OMEGA = 1.            *      
!  X      : N-vector X(1:N) that contains the solution vector    *      
!  RES    : N-vector RES(1:N) containing the residuum B - AX;    *      
!           the residuum is available even if the desired        *      
!           accuracy EPS could not be achieved with the given    *      
!           maximum number of iterations.                        *      
!  ITNUMB : num,bert of iterations actually performed            *      
!  IERR   : error parameter:                                     *      
!           = 0, the desired convergence criterium has not been  *      
!                met                                             *      
!           = 1, the solution X has been found                   *      
!           = 2, the desired accuracy has not been achieved after*      
!                KMAX iterations                                 *      
!           = 3, input data incorrect                            *      
!           = 4, system matrix A is numerically singular         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: GAUSEI, MNORM, CONV, RESID, MACHPD      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Gisela Engeln-MÅllges                             *      
!  Date      : 06.09.1992                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  Declarations                                                         
!                                                                       
      DOUBLEPRECISION A (1:IA, 1:N), B (1:N), X (1:N), WORK (1:N, 1:3), &
      RES (1:N), EPS, OMEGA, FMACHP, HELP, DIFFN, Q, RELERR, SUM, XN    
!                                                                       
!  Checking the inputs EPS, KMAX, IMETH and ISWITC                      
!                                                                       
      IF (EPS.LE.0.0D0.OR.KMAX.LT.1.OR.ISWITC.LT.0.OR.ISWITC.GT.3.OR.IME&
     &TH.LT.0.OR.IMETH.GT.2) THEN                                       
         IERR = 3 
         RETURN 
      ENDIF 
!                                                                       
!  Initialize the parameters KADAPT and OMEGA depending on the method   
!                                                                       
      IF (IMETH.EQ.0) THEN 
         OMEGA = 1.0D0 
      ELSEIF (IMETH.EQ.1) THEN 
         KADAPT = KMAX 
      ELSEIF (IMETH.EQ.2) THEN 
         KADAPT = KMAX 
         OMEGA = 1.0D0 
      ENDIF 
!                                                                       
!  Compute the machine constant and initialize the relative error bound 
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      RELERR = FMACHP * 8.0D0 
!                                                                       
!  Initialize                                                           
!                                                                       
      Q = 1.0D0 
      ITNUMB = 0 
!                                                                       
!  Check whether A is singular; if so, set IERR = 4.                    
!                                                                       
      DO 20 I = 1, N 
         SUM = DABS (A (I, 1) ) 
         DO 30 K = 2, N 
            SUM = SUM + DABS (A (I, K) ) 
   30    END DO 
         IF (SUM.EQ.0.0D0) THEN 
            IERR = 4 
            RETURN 
         ELSEIF (DABS (A (I, I) ) / SUM.LT.RELERR) THEN 
            IERR = 4 
            RETURN 
         ENDIF 
   20 END DO 
!                                                                       
!  Redefine the entries in A and B:   A(I,J) := A(I,J)/A(I,I)           
!  and B(I) := B(I)/A(I,I) .                                            
!                                                                       
      DO 40 I = 1, N 
         HELP = 1.0D0 / A (I, I) 
         DO 50 J = 1, N 
            A (I, J) = A (I, J) * HELP 
   50    END DO 
         B (I) = B (I) * HELP 
   40 END DO 
!                                                                       
!  Check for convergence                                                
!                                                                       
      IF (ISWITC.NE.0) THEN 
         CALL CONV (ISWITC, A, N, IA, IERR) 
         IF (IERR.EQ.0) RETURN 
      ENDIF 
!                                                                       
!  The vector RES serves as auxiliary storage for the previous solution 
!  vektor. Initially RES contains the staring vector.                   
!                                                                       
      DO 60 I = 1, N 
         RES (I) = X (I) 
   60 END DO 
!                                                                       
!  One iteration with the Gau·-Seidel method gives the first iterate X  
!                                                                       
      CALL GAUSEI (A, N, IA, B, OMEGA, X) 
!                                                                       
!  Up the iteration counter                                             
!                                                                       
      ITNUMB = ITNUMB + 1 
!                                                                       
!  Compute the difference of the last two iterates                      
!                                                                       
      DO 70 I = 1, N 
         WORK (I, 1) = X (I) - RES (I) 
   70 END DO 
!                                                                       
!  Iteration loop for the chosen method                                 
!                                                                       
      DO 80 K = 1, KMAX - 1 
!                                                                       
!  Check break-off criterion                                            
!                                                                       
         CALL MNORM (WORK (1, 1), N, DIFFN) 
         CALL MNORM (X, N, XN) 
         IF (DIFFN.LE.EPS * XN) THEN 
            IERR = 1 
            ITNUMB = K 
            CALL RESID (A, N, IA, B, X, RES) 
            RETURN 
         ENDIF 
         IF (K.EQ.KMAX - 1) THEN 
            ITNUMB = KMAX 
            IERR = 2 
            CALL RESID (A, N, IA, B, X, RES) 
            RETURN 
         ENDIF 
!                                                                       
!  RES contains the previous iterate                                    
!                                                                       
         DO 90 I = 1, N 
            RES (I) = X (I) 
   90    END DO 
!                                                                       
!  One iteration step using Gau·-Seidel for a fixed OMEGA               
!                                                                       
         CALL GAUSEI (A, N, IA, B, OMEGA, X) 
!                                                                       
!  Compute the difference of the last two iterates                      
!                                                                       
         DO 100 I = 1, N 
            WORK (I, 2) = X (I) - RES (I) 
  100    END DO 
!                                                                       
!  If the number of performed iterations K is divisible by KADAPT,      
!  then we compute Q in order to adjust the relaxation parameter;       
!  Q is an estimate of the spectral radius of the iteration matrix.     
!                                                                       
         IF (MOD (K, KADAPT) .EQ.0) THEN 
            DO 110 I = 1, N 
               IF (DABS (WORK (I, 1) ) .LT.FMACHP) THEN 
                  WORK (I, 3) = 1.0D0 
               ELSE 
                  WORK (I, 3) = WORK (I, 2) / WORK (I, 1) 
               ENDIF 
  110       END DO 
            CALL MNORM (WORK (1, 3), N, Q) 
!                                                                       
!  If Q > 1, then the iteration counter is upped by one and             
!  the next Gau·-Seidel step is executed; otherwise a new               
!  relaxation parameter is calculated.                                  
!                                                                       
            IF (Q.LE.1.0D0) THEN 
               Q = MAX (Q, OMEGA - 1.0D0) 
               OMEGA = 2.0D0 / (1.0D0 + DSQRT (1.0D0 - ( (Q + OMEGA -   &
               1.0D0) / OMEGA) **2 / Q) )                               
            ENDIF 
         ENDIF 
!                                                                       
!  The difference vector of the last two iterations is replaced         
!  by the one of the previous two iterations for the approximate solutio
!                                                                       
         DO 120 I = 1, N 
            WORK (I, 1) = WORK (I, 2) 
  120    END DO 
   80 END DO 
      END SUBROUTINE ADSOR                          
!                                                                       
!                                                                       
      SUBROUTINE GAUSEI (A, N, IA, B, OMEGA, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine performs one iteration with the Gau·-Seidel   *      
!  method for a given relaxation parameter.                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : 2-dim. array A(1:IA, 1:N), that contains the         *      
!           modified system matrix A : A(I,J)=A(I,J)/A(I,I) for  *      
!           I,J=1, ..., N                                        *      
!  N      : order of the system                                  *      
!  IA     : leading dimension of A, as specified in the calling  *      
!           program                                              *      
!  B      : N-vector B(1:N) with the modified right hand side:   *      
!           B(I)=B(I)/A(I,I); I=1, ..., N                        *      
!  OMEGA  : relaxation parameter                                 *      
!  X      : N-vector X(1:N) containing the starting vector for   *      
!           the iteration                                        *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X      : N-vector X(1:N) containing the next iteration vector *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gisela Engeln-MÅllges                              *      
!  Date     : 06.09.1992                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION A (1:IA, 1:N), B (1:N), X (1:N), OMEGA, S 
!                                                                       
      DO 10 I = 1, N 
         S = B (I) 
         DO 20 J = 1, N 
            S = S - A (I, J) * X (J) 
   20    END DO 
         X (I) = X (I) + OMEGA * S 
   10 END DO 
      RETURN 
      END SUBROUTINE GAUSEI                         
!                                                                       
!                                                                       
      SUBROUTINE MNORM (X, N, XNORM) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine calculates the maximum norm XNORM of an       *      
!  N-vector X.                                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gisela Engeln-MÅllges                              *      
!  Date     : 06.09.1992                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X (1:N), XNORM 
!                                                                       
      XNORM = DABS (X (1) ) 
      DO 10 I = 2, N 
         XNORM = DMAX1 (XNORM, DABS (X (I) ) ) 
   10 END DO 
      RETURN 
      END SUBROUTINE MNORM                          
!                                                                       
!                                                                       
      SUBROUTINE CONV (ISWITC, A, N, IA, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine helps check convergence.                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  ISWITC : Parameter that determines the convergence criterion  *      
!           to be checked:                                       *      
!           = 0, none                                            *      
!           = 1, row sum criterion                               *      
!           = 2, column sum criterion                            *      
!           = 3, criterion of Schmidt and v. Mises               *      
!  A      : 2-dim. array A(1:IA, 1:N), containing the matrix for *      
!           which we want to check convergence of the iterates   *      
!           from the various SOR algorithms                      *      
!  N      : order of the matrix A                                *      
!  IA     : leading dimension of A, as prescribed in the calling *      
!           program                                              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  IERR   : error parameter:                                     *      
!           = 0, the desired convergence criterion has not been  *      
!                met                                             *      
!           = 1, the desired criterion is satified               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gisela Engeln-MÅllges                              *      
!  Date     : 06.09.1992                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION A (1:IA, 1:N), SUM 
!                                                                       
!  Row sum criterion                                                    
!                                                                       
      IF (ISWITC.EQ.1) THEN 
         DO 10 I = 1, N 
            SUM = - 1.0D0 
            DO 20 J = 1, N 
               SUM = SUM + DABS (A (I, J) ) 
   20       END DO 
            IF (SUM.LT.1.0D0) THEN 
               IERR = 1 
            ELSE 
               IERR = 0 
               RETURN 
            ENDIF 
   10    END DO 
!                                                                       
!  Column sum criterion                                                 
!                                                                       
      ELSEIF (ISWITC.EQ.2) THEN 
         DO 30 J = 1, N 
            SUM = - 1.0D0 
            DO 40 I = 1, N 
               SUM = SUM + DABS (A (I, J) ) 
   40       END DO 
            IF (SUM.LT.1.0D0) THEN 
               IERR = 1 
            ELSE 
               IERR = 0 
               RETURN 
            ENDIF 
   30    END DO 
!                                                                       
!  Criterion of Schmidt and v. Mises                                    
!                                                                       
      ELSEIF (ISWITC.EQ.3) THEN 
         SUM = - N 
         DO 50 I = 1, N 
            DO 60 J = 1, N 
               SUM = SUM + A (I, J) * A (I, J) 
   60       END DO 
   50    END DO 
         SUM = DSQRT (SUM) 
         IF (SUM.LT.1.0D0) THEN 
            IERR = 1 
         ELSE 
            IERR = 0 
            RETURN 
         ENDIF 
      ENDIF 
      END SUBROUTINE CONV                           
!                                                                       
!                                                                       
      SUBROUTINE RESID (A, N, IA, B, X, RES) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine computes the residuum  RES = B - AX, where    *      
!  both A and B are given in modified form.                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gisela Engeln-MÅllges                              *      
!  Date     : 09.06.1992                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION A (1:IA, 1:N), B (1:N), X (1:N), RES (1:N),       &
      DSUM                                                              
!                                                                       
      DO 10 I = 1, N 
         DSUM = B (I) 
         DO 20 J = 1, N 
            DSUM = DSUM - A (I, J) * X (J) 
   20    END DO 
         RES (I) = DSUM 
   10 END DO 
      RETURN 
      END SUBROUTINE RESID                          

![  {Systems with Five--Diagonal Symmetric Matrices}                    
![  {Systems with Five--Diagonal Symmetric Matrices}*)                  
      SUBROUTINE FDISY (N, DM, DU1, DU2, RS, X, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Solving a system of linear equations                          *      
!                     A * X = RS                                 *      
!  for a five-diagonal, symmetric and strongly nonsingular       *      
!  matrix A. The matrix A is given by the three N-vectors DM,    *      
!  DU1 and DU2. The system of equations has the form :           *      
!                                                                *      
!  DM(1)*X(1) + DU1(1)*X(2) + DU2(1)*X(3)               = RS(1)  *      
!  DU1(1)*X(1) + DM(2)*X(2) + DU1(2)*X(3) + DU2(2)*X(4) = RS(2)  *      
!                                                                *      
!  DU2(I-2)*X(I-2) + DU1(I-1)*X(I-1) + DM(I)*X(I) +              *      
!                       + DU1(I)*X(I+1) + DU2(I)*X(I+2) = RS(I)  *      
!             for I = 3, ..., N - 2, and                         *      
!                                                                *      
!  DU2(N-3)*X(N-2) + DU1(N-2)*X(N-1) + DM(N-1)*X(N-1) +          *      
!                                       + DU1(N-1)*X(N) = RS(N-1)*      
!  DU2(N-2)*X(N-2) + OD(N-1)*X(N-1) + DM(N)*X(N)        = RS(N)  *      
!                                                                *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    : number of equations, N > 3                             *      
!  DM   : N-vector DM(1:N); main diagonal of A                   *      
!         DM(1), DM(2), ... , DM(N)                              *      
!  DU1  : N-vector DU1(1:N); co-diagonal of A                    *      
!         DU1(1), DU1(2), ... , DU1(N-1)                         *      
!  DU2  : N-vector DU2(1:N); second co-diagonal of A             *      
!         DU2(1), DU2(2), ... , DU2(N-2)                         *      
!  RS   : N-vector RS(1:N); the right hand side                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  DM   :)                                                       *      
!  DU1  :) overwritten with intermediate quantities              *      
!  DU2  :)                                                       *      
!  RS   :)                                                       *      
!  X    : N-vector X(1:N) containing the solution vector         *      
!  MARK : error parameter                                        *      
!         MARK=-2 : condition N > 3 is not satisfied             *      
!         MARK=-1 : A is strongly nonsingular, but not positive  *      
!                   definite                                     *      
!         MARK= 0 : numerically the matrix A is not strongly     *      
!                   nonsingular                                  *      
!         MARK= 1 : A is positive definite                       *      
!                                                                *      
!  NOTE: If MARK = +/- 1, then the determinant of A is:          *      
!           DET A = DM(1) * DM(2) * ... * DM(N)                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: FDISYP, FDISYS, MACHPD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors  : Gisela Engeln-Muellges                             *      
!  date     : 01.07.1992                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (1:N), DU1 (1:N), DU2 (1:N), RS (1:N), X (1:N) 
      MARK = - 2 
      IF (N.LT.4) RETURN 
!                                                                       
!  Factorization of the matrix A                                        
!                                                                       
      CALL FDISYP (N, DM, DU1, DU2, MARK) 
!                                                                       
!  if MARK = +/- 1 , update and backsubstitute                          
!                                                                       
      IF (MARK.EQ.1) THEN 
         CALL FDISYS (N, DM, DU1, DU2, RS, X) 
      ENDIF 
      RETURN 
      END SUBROUTINE FDISY                          
!                                                                       
!                                                                       
      SUBROUTINE FDISYP (N, DM, DU1, DU2, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Factor a five-diagonal, symmetric and strongly nonsingular    *      
!  matrix A, that is given by the three N-vectors DM, DU1 and    *      
!  DU2, into its Cholesky factors A =  R(TRANSP) * D * R  by     *      
!  applying the root-free Cholesky method for five-diagonal      *      
!  matrices. The form of the linear system is identical with     *      
!  the one in SUBROUTINE FDISY.                                  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    : number of equations, N > 3                             *      
!  DM   : N-vector DM(1:N); main diagonal of A                   *      
!         DM(1), DM(2), ... , DM(N)                              *      
!  DU1  : N-vector DU1(1:N); upper co-diagonal of A              *      
!         DU1(1), DU1(2), ... , DU1(N-1)                         *      
!  DU2  : N-vector DU2(1:N); second upper co-diagonal of A       *      
!         DU2(1), DU2(2), ... , DU2(N-2);                        *      
!         due to symmetry the lower co-diagonals do not need to  *      
!         be stored separately.                                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  DM   :) overwritten with auxiliary vectors containing the     *      
!  DU1  :) Cholesky factors of A. The co-diagonals of the unit   *      
!  DU2  :) upper tridiagonal matrix R are stored in DU1 and DU2, *      
!          the diagonal matrix D in DM.                          *      
!  MARK : error parameter                                        *      
!         MARK=-2 : condition N > 3 is not satisfied             *      
!         MARK=-1 : A is strongly nonsingular, but not positive  *      
!                   definite                                     *      
!         MARK= 0 : numerically the matrix is not strongly       *      
!                   nonsingular                                  *      
!         MARK= 1 : A is positive definite                       *      
!                                                                *      
!  NOTE : If MARK = +/-1, then the inertia of A, i. e., the      *      
!         number of positive and negative eigenvalues of A,      *      
!         is the same as the number of positive and negative     *      
!         numbers among the components of DM.                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors  : Gisela Engeln-Muellges                             *      
!  date     : 01.07.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (1:N), DU1 (1:N), DU2 (1:N) 
!                                                                       
!   calculating the machine constant                                    
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!   determining the relative error bound                                
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!   checking for N > 3                                                  
!                                                                       
      MARK = - 2 
      IF (N.LT.4) RETURN 
      DU1 (N) = 0.0D0 
      DU2 (N) = 0.0D0 
      DU2 (N - 1) = 0.0D0 
!                                                                       
!   checking for strong nonsingularity of the matrix A for N=1          
!                                                                       
      ROW = DABS (DM (1) ) + DABS (DU1 (1) ) + DABS (DU2 (1) ) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      IF (DM (1) .LT.0.0D0) THEN 
         MARK = - 1 
         RETURN 
      ELSEIF (DABS (DM (1) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!                                                                       
!   factoring A while checking for strong nonsingularity                
!                                                                       
      DUMMY = DU1 (1) 
      DU1 (1) = DU1 (1) / DM (1) 
      DUMMY1 = DU2 (1) 
      DU2 (1) = DU2 (1) / DM (1) 
      ROW = DABS (DUMMY) + DABS (DM (2) ) + DABS (DU1 (2) ) + DABS (DU2 &
      (2) )                                                             
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      DM (2) = DM (2) - DUMMY * DU1 (1) 
      IF (DM (2) .LT.0.0D0) THEN 
         MARK = - 1 
         RETURN 
      ELSEIF (DABS (DM (2) ) .LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      DUMMY = DU1 (2) 
      DU1 (2) = (DU1 (2) - DUMMY1 * DU1 (1) ) / DM (2) 
      DUMMY2 = DU2 (2) 
      DU2 (2) = DU2 (2) / DM (2) 
      DO 20 I = 3, N, 1 
         ROW = DABS (DUMMY1) + DABS (DUMMY) + DABS (DM (I) ) + DABS (   &
         DU1 (I) ) + DABS (DU2 (I) )                                    
         IF (ROW.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         D = 1.0D0 / ROW 
         DM (I) = DM (I) - DM (I - 1) * DU1 (I - 1) * DU1 (I - 1)       &
         - DUMMY1 * DU2 (I - 2)                                         
         IF (DM (I) .LT.0.0D0) THEN 
            MARK = - 1 
            RETURN 
         ELSEIF (DABS (DM (I) ) * D.LE.EPS) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         IF (I.LT.N) THEN 
            DUMMY = DU1 (I) 
            DU1 (I) = (DU1 (I) - DUMMY2 * DU1 (I - 1) ) / DM (I) 
            DUMMY1 = DUMMY2 
         ENDIF 
         IF (I.LT.N - 1) THEN 
            DUMMY2 = DU2 (I) 
            DU2 (I) = DU2 (I) / DM (I) 
         ENDIF 
   20 END DO 
      MARK = 1 
      RETURN 
      END SUBROUTINE FDISYP                         
!                                                                       
!                                                                       
      SUBROUTINE FDISYS (N, DM, DU1, DU2, RS, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Solving a linear system of equations                          *      
!               A * X = RS                                       *      
!  for a five-diagonal, symmetric and strongly nonsingular       *      
!  matrix A, once its Cholesky factors have been calculated by   *      
!  SUBROUTINE FDISYP. Here the factors of A are used as input    *      
!  arrays and they are stored in the three N-vectors DM, DU1     *      
!  and DU2.                                                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETER:                                              *      
!  ================                                              *      
!  N    : number of equations, N > 3                             *      
!  DM   : N-vector DM(1:N);  diagonal matrix D                   *      
!  DU1  : N-vector DM(1:N); ) co-diagonals of the upper          *      
!  DU2  : N-vector DM(1:N); ) triangular  matrix R               *      
!  RS   : N-vector DM(1:N); the right hand side                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  X    : N-vector X(1:N) containing the solution vector         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 29.04.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (1:N), DU1 (1:N), DU2 (1:N), RS (1:N), X (1:N) 
!                                                                       
!  updating                                                             
!                                                                       
      DUMMY1 = RS (1) 
      RS (1) = DUMMY1 / DM (1) 
      DUMMY2 = RS (2) - DU1 (1) * DUMMY1 
      RS (2) = DUMMY2 / DM (2) 
      DO 10 I = 3, N, 1 
         DUMMY1 = RS (I) - DU1 (I - 1) * DUMMY2 - DU2 (I - 2) * DUMMY1 
         RS (I) = DUMMY1 / DM (I) 
         DUMMY3 = DUMMY2 
         DUMMY2 = DUMMY1 
         DUMMY1 = DUMMY3 
   10 END DO 
!                                                                       
!  backsubstitution                                                     
!                                                                       
      X (N) = RS (N) 
      X (N - 1) = RS (N - 1) - DU1 (N - 1) * X (N) 
      DO 20 I = N - 2, 1, - 1 
         X (I) = RS (I) - DU1 (I) * X (I + 1) - DU2 (I) * X (I + 2) 
   20 END DO 
      RETURN 
      END SUBROUTINE FDISYS                         

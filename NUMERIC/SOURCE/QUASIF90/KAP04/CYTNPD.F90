      SUBROUTINE CYTNPD (N, DM, DU, CR, RS, X, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a system of linear  equations                      *      
!                  A * X = RS                                    *      
!     for a cyclically tridiagonal, symmetric, strongly          *      
!     nonsingular matrix A. The matrix A is given by two         *      
!     N-vectors DM and DU. The system of equations has the form: *      
!                                                                *      
!     DM(1) * X(1) + DU(1) * X(2) + DU(N) * X(N)       = RS(1)   *      
!                                                                *      
!     DU(I-1) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1) = RS(I)   *      
!                for I = 2, ..., N - 1, and                      *      
!                                                                *      
!     DU(N) * X(1) + DU(N-1) * X(N-1) + DM(N) * X(N)   = RS(N)   *      
!                                                                *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DM   : N-vector DM(1:N); main diagonal of A                *      
!            DM(1), DM(2), ... , DM(N)                           *      
!     DU   : N-vector DU(1:N); upper co-diagonal of A            *      
!            DU(1), DU(2), ... , DU(N-1); the off-diagonal       *      
!            element A(1,N) is stored in DU(N).                  *      
!     RS   : N-vector RS(1:N); the right hand side               *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DM   :)                                                    *      
!     DU   :) overwritten with intermediate vectors              *      
!     CR   :)                                                    *      
!     RS   :)                                                    *      
!     X    : N-vector X(1:N), containing the solution            *      
!     MARK : error parameter                                     *      
!            MARK=-2 : condition N > 2 is not satisfied          *      
!            MARK=-1 : A is strongly nonsingular, but not        *      
!                      positive definite                         *      
!            MARK= 0 : numerically the matrix A is not strongly  *      
!                      nonsingular                               *      
!            MARK= 1 : A is positive definite                    *      
!                                                                *      
!     NOTE : If MARK = +/- 1, the determinant of A is given as:  *      
!               DET A = DM(1) * DM(2) * ... * DM(N)              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: CYTSYP, CYTSYS, MACHPD                  *      
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
      DOUBLEPRECISION DM (1:N), DU (1:N), CR (1:N), RS (1:N), X (1:N) 
      MARK = - 2 
      IF (N.LT.3) RETURN 
!                                                                       
!  factorization of the matrix A                                        
!                                                                       
      CALL CNPSYP (N, DM, DU, CR, MARK) 
!                                                                       
!  if MARK = +/- 1, update and backsubstitute                           
!                                                                       
      IF ( (MARK.EQ.1) .OR. (MARK.EQ. - 1) ) THEN 
         CALL CYTSYS (N, DM, DU, CR, RS, X) 
      ENDIF 
      RETURN 
      END SUBROUTINE CYTNPD                         
!                                                                       
!                                                                       
      SUBROUTINE CNPSYP (N, DM, DU, CR, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factoring a cyclically tridiagonal, symmetric and strongly *      
!     nonsingular matrix A, that is given by the two N-vectors   *      
!     DM and DU, into its Cholesky factors                       *      
!                    A = R(TRANSP) * D * R                       *      
!     by applying the root-free Cholesky-method for tridiagonal  *      
!     cyclic matrices. The form of the system of equations is    *      
!     identical to the one described in SUBROUTINE CYTSY.        *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DM   : N-vector DM(1:N); main diagonal of A                *      
!            DM(1), DM(2), ... , DM(N)                           *      
!     DU   : N-vector DU(1:N); upper co-diagonal of A            *      
!            DU(1), DU(2), ... , DU(N-1); the off-diagonal       *      
!            element A(1,N) is stored in DU(N).                  *      
!            Due to symmetry the lower co-diagonal does not need *      
!            to be stored separately.                            *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     DM   :) overwritten with auxiliary vectors from the        *      
!     DU   :) factorization of A. The co-diagonal of the unit    *      
!     CR   :) upper tridiagonal matrix R is stored in DU, the    *      
!             diagonal matrix D appears in DM and the right hand *      
!             side in CR.                                        *      
!     MARK : error parameter                                     *      
!            MARK=-2 : condition N > 2 is not satisfied          *      
!            MARK=-1 : A is strongly nonsingular, but not        *      
!                      positive definite                         *      
!            MARK= 0 : numerically A is not strongly             *      
!                      nonsingular                               *      
!            MARK= 1 : A is positive definite                    *      
!                                                                *      
!     NOTE : If MARK = +/- 1, then the inertia of A, i. e., the  *      
!            number of positive and negative eigenvalues of A,   *      
!            is the same as the number of positive and negative  *      
!            numbers among the components of DM.                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutine required: MACHPD                                   *      
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
      DOUBLEPRECISION DM (1:N), DU (1:N), CR (1:N) 
!                                                                       
!   calculating the machine constant                                    
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!   determinaing the relative error bound                               
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!   testing of condition N > 2                                          
!                                                                       
      MARK = - 2 
      IF (N.LT.3) RETURN 
      MARK = 1 
!                                                                       
!   checking for strong nonsingularity of A for N=1                     
!                                                                       
      ROW = DABS (DM (1) ) + DABS (DU (1) ) + DABS (DU (N) ) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      IF (DM (1) .LT.0.0D0) MARK = - 1 
      IF (DABS (DM (1) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!                                                                       
!   factoring A while checking for strong nonsingularity                
!                                                                       
      DUMMY = DU (1) 
      DU (1) = DU (1) / DM (1) 
      CR (1) = DU (N) / DM (1) 
      DO 20 I = 2, N - 1, 1 
         ROW = DABS (DM (I) ) + DABS (DU (I) ) + DABS (DUMMY) 
         IF (ROW.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         D = 1.0D0 / ROW 
         DM (I) = DM (I) - DUMMY * DU (I - 1) 
         IF (DM (I) .LT.0.0D0) MARK = - 1 
         IF (DABS (DM (I) ) * D.LE.EPS) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         IF (I.LT. (N - 1) ) THEN 
            CR (I) = - DUMMY * CR (I - 1) / DM (I) 
            DUMMY = DU (I) 
            DU (I) = DU (I) / DM (I) 
         ELSE 
            DUMMY2 = DU (I) 
            DU (I) = (DU (I) - DUMMY * CR (I - 1) ) / DM (I) 
         ENDIF 
   20 END DO 
      ROW = DABS (DU (N) ) + DABS (DM (N) ) + DABS (DUMMY2) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      DM (N) = DM (N) - DM (N - 1) * DU (N - 1) * DU (N - 1) 
      DUMMY = 0.0D0 
      DO 30 I = 1, N - 2, 1 
         DUMMY = DUMMY + DM (I) * CR (I) * CR (I) 
   30 END DO 
      DM (N) = DM (N) - DUMMY 
      IF (DM (N) .LT.0) MARK = - 1 
      IF (DABS (DM (N) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      RETURN 
      END SUBROUTINE CNPSYP                         

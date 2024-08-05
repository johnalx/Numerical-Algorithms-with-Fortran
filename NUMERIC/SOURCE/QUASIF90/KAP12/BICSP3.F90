      SUBROUTINE BICSP3 (N, M, A, X, Y, FN, F, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Determining bicubic splines for given functional values and   *      
!  surface normals at all points.                                *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    : number of X-intervals                                  *      
!  M    : number of Y-intervals                                  *      
!  A    : 4-dimensional array A(0:N,0:M,0:KDIM,0:LDIM); contains *      
!         the spline coefficients. On call, A(I,J,0,0) must      *      
!         contain the functional values U(I,J).                  *      
!         BICSP1 determines all other A(I,J,K,L) for I=0 to N-1  *      
!         and J=0 to M-1. Elements A(N,M,K,L), that are not      *      
!         assigned a value on call, remain unassigned.           *      
!  X    : (N+1)-vector X(0:N) containing the endpoints of the    *      
!         X-intervals                                            *      
!  Y    : (N+1)-vector Y(0:M) containing the endpoints of the    *      
!         Y-intervals                                            *      
!  FN   : 3-dimensional array FN(0:N,0:M,1:3) containing the     *      
!         normal vectors at all points.                          *      
!  F    : auxiliary vector F(1:9*MAX(N,M)-8)                     *      
!  IERR : is initially set to 0. Will be set different from zero *      
!         if the algorithm detects an error. If errors occur     *      
!         the program does not complete the calculations.        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BIC3S1, BIC2S3, BIC1S3, BIC1S4, BIC1S5, *      
!                        BIC1S6, BIC1S7, BIC1S8, BIC1S9, TRIDIG  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), Y (0:M), F ( * ) 
      DIMENSION AA (0:3, 0:3), FN (0:N, 0:M, 1:3) 
!                                                                       
!*  steps 1 and 2                                                       
!                                                                       
      CALL BIC3S1 (N, M, A, FN, JERR) 
      IF (JERR.NE.0) RETURN 
      IERR = JERR 
!                                                                       
!*  splitting of vector F for finding X values                          
!                                                                       
      IH = 1 
      ISA = IH + N 
      ISB = ISA + N - 1 
      ISC = ISB + N - 1 
      ISD = ISC + N - 1 
      ISX = ISD+N - 1 
      ISGAMM = ISX + N - 1 
      ISALPH = ISGAMM + N - 1 
      ISG = ISALPH + N - 1 
!                                                                       
!*  splitting of the vector F for finding Y values                      
!                                                                       
      JH = 1 
      JSA = JH + M 
      JSB = JSA + M - 1 
      JSC = JSB + M - 1 
      JSD = JSC + M - 1 
      JSX = JSD+M - 1 
      JSGAMM = JSX + M - 1 
      JSALPH = JSGAMM + M - 1 
      JSG = JSALPH + M - 1 
!                                                                       
!*  step 3                                                              
!                                                                       
      CALL BIC2S3 (N, M, A, X, F (IH), JERR) 
      IERR = JERR + 4 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 3 continued                                                    
!                                                                       
      CALL BIC1S3 (N, M, A, X, F (IH), F (ISA), F (ISB), F (ISC),       &
      F (ISD), F (ISX), F (ISGAMM), F (ISALPH), F (ISG), JERR)          
      IERR = JERR + 4 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 4                                                              
!                                                                       
      CALL BIC1S4 (N, M, A, Y, F (JH), F (JSA), F (JSB), F (JSC),       &
      F (JSD), F (JSX), F (JSGAMM), F (JSALPH), F (JSG), JERR)          
      IERR = JERR + 6 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  steps 5, 6, 7 are contained in step 8                               
!*  loops over all X and Y values                                       
!                                                                       
      DO 112 I = 0, N - 1 
         DO 111 J = 0, M - 1 
            CALL BIC1S8 (N, M, A, X, Y, I, J, AA) 
!                                                                       
!*  transfer  AA to A                                                   
!                                                                       
            CALL BIC1S9 (N, M, A, I, J, AA) 
  111    END DO 
  112 END DO 
!                                                                       
!*  all spline-coefficients are determined now                          
!                                                                       
      IERR = 0 
      RETURN 
      END SUBROUTINE BICSP3                         
!                                                                       
!                                                                       
      SUBROUTINE BIC3S1 (N, M, A, FN, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  step 1 and step 2:                                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), FN (0:N, 0:M, 3) 
!                                                                       
!*  The array FN contains of the normal vectors                         
!                                                                       
      DO 102 I = 0, N 
         DO 101 J = 0, M 
            IF (FN (I, J, 3) .EQ.0.0D0) GOTO 900 
            A (I, J, 1, 0) = - FN (I, J, 1) / FN (I, J, 3) 
            A (I, J, 0, 1) = - FN (I, J, 2) / FN (I, J, 3) 
  101    END DO 
  102 END DO 
      IERR = 0 
      RETURN 
!                                                                       
!*  error, third component of one normal vector is zero                 
!                                                                       
  900 IERR = 1 
      RETURN 
      END SUBROUTINE BIC3S1                         

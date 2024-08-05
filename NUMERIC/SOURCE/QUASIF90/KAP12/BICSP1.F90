![KA{P 12}{Two--Dim., B\'ezier, Surface, B--Splines}                    
![        {Two--Dimensional Splines, Surface Splines,                   
![         B\'ezier Splines, B--Splines}*)                              
![  {Interpolating Two--Dimensional Cubic Splines}                      
![  {Interpolating Two--Dimensional Cubic Splines for                   
![   Constructing Smooth Surfaces}*)                                    
      SUBROUTINE BICSP1 (N, M, A, X, Y, F, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Determination of bicubic splines.                             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    : number of X-intervals                                  *      
!  M    : number of Y-intervals                                  *      
!  A    : 4-dimensional array A(0:N,0:M,0:KDIM,0:LDIM); contains *      
!                     the spline coefficients: on call,          *      
!          A(I,J,0,0) must contain the functional values U(I,J)  *      
!          A(I,J,1,0) must contain the derivatives P(I,J) for    *      
!                     J=0 and J=M                                *      
!          A(I,J,0,1) must contain the derivatives Q(I,J) for    *      
!                     I=0 and I=N                                *      
!          A(I,J,1,1) must contain the derivatives R(I,J) for    *      
!                     I=0 and I=N, as well as for J=0 and J=M    *      
!         BICSP1 determines all other A(I,J,K,L) for I=0 to N-1  *      
!         and J=0 to M-1. Elements A(N,M,K,L), that are unas-    *      
!         signed on call, will remain unassigned.                *      
!  X    : (N+1)-vector X(0:N) containing the endpoints of the    *      
!         X-intervals                                            *      
!  Y    : (N+1)-vector Y(0:M) containing the endpoints of the    *      
!         Y-intervals                                            *      
!  F    : auxiliary vector F(1:9*MAX(N,M)-8)                     *      
!  IERR : initially set to 0. Will be set different from zero,   *      
!         if the algorithm detects an error. If errors occur the *      
!         program does not complete the computations.            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BIC1S1, BIC1S2, BIC1S3, BIC1S4, BIC1S5, *      
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
      DIMENSION AA (0:3, 0:3) 
!                                                                       
!*  splitting the vector F for finding X values                         
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
!*  splitting the vector F for finding Y values                         
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
!*  step 1                                                              
!                                                                       
      CALL BIC1S1 (N, M, A, X, F (IH), F (ISA), F (ISB), F (ISC),       &
      F (ISD), F (ISX), F (ISGAMM), F (ISALPH), F (ISG), JERR)          
      IERR = JERR 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 2                                                              
!                                                                       
      CALL BIC1S2 (N, M, A, Y, F (JH), F (JSA), F (JSB), F (JSC),       &
      F (JSD), F (JSX), F (JSGAMM), F (JSALPH), F (JSG), JERR)          
      IERR = JERR + 2 
      IF (JERR.NE.0) RETURN 
!                                                                       
!*  step 3                                                              
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
!*  steps 5, 6, 7 are contained in step 8 ;                             
!*  looping over all X-Y values                                         
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
!*  all spline-coefficients have now been determined                    
!                                                                       
      IERR = 0 
      RETURN 
      END SUBROUTINE BICSP1                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S1 (N, M, A, X, H, SA, SB, SC, SD, SX, SGAMM,      &
      SALPH, SG, IERR)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
!  step 1:                                                       *      
!                                                                *      
!  All vectors, that are needed for solving the linear system of *      
!  equations start with the letter S. The notation is taken      *      
!  directly from the formulas of the text section, chapter 12.1. *      
!  For a description of the parameters see BICSP1.               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRIDIG                                  *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N) 
      DIMENSION H (0:N - 1) 
      DIMENSION SX (N - 1), SA (N - 1), SB (N - 1), SC (N - 1), SG (N - &
      1), SGAMM (N - 1)                                                 
      DIMENSION SALPH (N - 1), SD (N - 1) 
!                                                                       
!*  determine H(I), test strict monotonicity of X(I)                    
!                                                                       
      DO 100 I = 0, N - 1 
         H (I) = X (I + 1) - X (I) 
         IF (H (I) .LE.0) GOTO 900 
  100 END DO 
!                                                                       
!*  determine the columns of the tridiagonal system                     
!                                                                       
      DO 101 I = 1, N - 1 
         SB (I) = 1.0D0 / H (I - 1) 
         SC (I) = 1.0D0 / H (I) 
         SD (I) = 2.0D0 * (SB (I) + SC (I) ) 
  101 END DO 
!                                                                       
!*  for the 1st system, TRIDIG is performed to the end                  
!                                                                       
      IREP = 0 
!                                                                       
!*  loop covering the M+1 systems of equations                          
!                                                                       
      DO 104 J = 0, M 
!                                                                       
!*  determine the right-hand side of each system                        
!                                                                       
         DO 102 I = 1, N - 1 
            SA (I) = 3.0D0 * ( (A (I, J, 0, 0) - A (I - 1, J, 0, 0) )   &
            / H (I - 1) **2 + (A (I + 1, J, 0, 0) - A (I, J, 0, 0) )    &
            / H (I) **2)                                                
  102    END DO 
!                                                                       
!*  corrections for the first and last equations only                   
!                                                                       
         SA (1) = SA (1) - A (0, J, 1, 0) / H (0) 
         SA (N - 1) = SA (N - 1) - A (N, J, 1, 0) / H (N - 1) 
!                                                                       
!*  solve the system, SX is the solution vector                         
!                                                                       
         CALL TRIDIG (N - 1, SA, SB, SC, SD, IREP, SX, IERR, SGAMM,     &
         SALPH, SG)                                                     
!                                                                       
!*  return with IERR=1, if the system is unsolvable                     
!                                                                       
         IF (IERR.NE.0) RETURN 
!                                                                       
!*  to increase speed from the 2nd system onwards (compare TRIDIG)      
!                                                                       
         IREP = 1 
!                                                                       
!*  transfer solution vector to array A                                 
!                                                                       
         DO 103 I = 1, N - 1 
            A (I, J, 1, 0) = SX (I) 
  103    END DO 
  104 END DO 
!                                                                       
!*  algorithm has run successfully                                      
!                                                                       
      RETURN 
!                                                                       
!*  return with IERR=2 if monotonicity of X(I) is not met               
!                                                                       
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC1S1                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S2 (N, M, A, Y, H, SA, SB, SC, SD, SX, SGAMM,      &
      SALPH, SG, IERR)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 2:                                                       *      
!                                                                *      
!  All vectors, that are needed for solving the system of        *      
!  equations start with the letter S. The notation is the same   *      
!  as that of the formulas in chapter 12.1.                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRIDIG                                  *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), Y (0:M) 
      DIMENSION H (0:M - 1) 
      DIMENSION SX (M - 1), SA (M - 1), SB (M - 1), SC (M - 1), SG (M - &
      1), SGAMM (M - 1)                                                 
      DIMENSION SALPH (M - 1), SD (M - 1) 
!                                                                       
!*  determine H(J), test strict monotonicity of Y(J)                    
!                                                                       
      DO 100 J = 0, M - 1 
         H (J) = Y (J + 1) - Y (J) 
         IF (H (J) .LE.0) GOTO 900 
  100 END DO 
!                                                                       
!*  determine the columns of the tridiagonal system                     
!                                                                       
      DO 101 J = 1, M - 1 
         SB (J) = 1.0D0 / H (J - 1) 
         SC (J) = 1.0D0 / H (J) 
         SD (J) = 2.0D0 * (SB (J) + SC (J) ) 
  101 END DO 
!                                                                       
!*  for the 1st system, TRIDIG is to be performed to the end            
!                                                                       
      IREP = 0 
!                                                                       
!*  loop covering the N+1 systems of equations                          
!                                                                       
      DO 104 I = 0, N 
!                                                                       
!*  determine the right-hand side of each system                        
!                                                                       
         DO 102 J = 1, M - 1 
            SA (J) = 3.0D0 * ( (A (I, J, 0, 0) - A (I, J - 1, 0, 0) )   &
            / H (J - 1) **2 + (A (I, J + 1, 0, 0) - A (I, J, 0, 0) )    &
            / H (J) **2)                                                
  102    END DO 
!                                                                       
!*  corrections in the first and last equations only                    
!                                                                       
         SA (1) = SA (1) - A (I, 0, 0, 1) / H (0) 
         SA (M - 1) = SA (M - 1) - A (I, M, 0, 1) / H (M - 1) 
!                                                                       
!*  solve the system, SX is the solution vector                         
!                                                                       
         CALL TRIDIG (M - 1, SA, SB, SC, SD, IREP, SX, IERR, SGAMM,     &
         SALPH, SG)                                                     
!                                                                       
!*  return with IERR=1, if the system is unsolvable                     
!                                                                       
         IF (IERR.NE.0) RETURN 
!                                                                       
!*  to increase the speed from the 2nd set onward (compare TRIDIG)      
!                                                                       
         IREP = 1 
!                                                                       
!*  transfer solution vector to array A                                 
!                                                                       
         DO 103 J = 1, M - 1 
            A (I, J, 0, 1) = SX (J) 
  103    END DO 
  104 END DO 
!                                                                       
!*  algorithm has run successfully                                      
!                                                                       
      RETURN 
!                                                                       
!*  return with IERR=2, if monotonicity of Y(J) is not true             
!                                                                       
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC1S2                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S3 (N, M, A, X, H, SA, SB, SC, SD, SX, SGAMM,      &
      SALPH, SG, IERR)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 3:                                                       *      
!                                                                *      
!  All vectors, that are used for solving the system of equations*      
!  start with the letter S. The notation follows directly from   *      
!  the formulas in chapter 12.1.                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRIDIG                                  *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N) 
      DIMENSION H (0:N - 1) 
      DIMENSION SX (N - 1), SA (N - 1), SB (N - 1), SC (N - 1), SG (N - &
      1), SGAMM (N - 1)                                                 
      DIMENSION SALPH (N - 1), SD (N - 1) 
!                                                                       
!*  determine H(I), test strict monotonicity of X(I)                    
!                                                                       
      DO 100 I = 0, N - 1 
         H (I) = X (I + 1) - X (I) 
         IF (H (I) .LE.0) GOTO 900 
  100 END DO 
!                                                                       
!*  determine the columns of the tridiagonal system                     
!                                                                       
      DO 101 I = 1, N - 1 
         SB (I) = 1.0D0 / H (I - 1) 
         SC (I) = 1.0D0 / H (I) 
         SD (I) = 2.0D0 * (SB (I) + SC (I) ) 
  101 END DO 
!                                                                       
!*  in the 1st system, TRIDIG is to be performed to the end             
!                                                                       
      IREP = 0 
!                                                                       
!*  loop covering all the systems of equations                          
!                                                                       
      DO 104 J = 0, M, M 
!                                                                       
!*  determine the right-hand side of each system                        
!                                                                       
         DO 102 I = 1, N - 1 
            SA (I) = 3.0D0 * ( (A (I, J, 0, 1) - A (I - 1, J, 0, 1) )   &
            / H (I - 1) **2 + (A (I + 1, J, 0, 1) - A (I, J, 0, 1) )    &
            / H (I) **2)                                                
  102    END DO 
!                                                                       
!*  corrections in the first and last equations only                    
!                                                                       
         SA (1) = SA (1) - A (0, J, 1, 1) / H (0) 
         SA (N - 1) = SA (N - 1) - A (N, J, 1, 1) / H (N - 1) 
!                                                                       
!*  solve the system, SX is the solution vector                         
!                                                                       
         CALL TRIDIG (N - 1, SA, SB, SC, SD, IREP, SX, IERR, SGAMM,     &
         SALPH, SG)                                                     
!                                                                       
!*  return with IERR=1, if the system is unsolvable                     
!                                                                       
         IF (IERR.NE.0) RETURN 
!                                                                       
!*  to increase the speed from the 2nd set onwards (compare TRIDIG)     
!                                                                       
         IREP = 1 
!                                                                       
!*  transfer solution vector to array A                                 
!                                                                       
         DO 103 I = 1, N - 1 
            A (I, J, 1, 1) = SX (I) 
  103    END DO 
  104 END DO 
!                                                                       
!*  algorithm has run successfully                                      
!                                                                       
      RETURN 
!                                                                       
!*  return with IERR=2, if monotonicity of X(I) is not true             
!                                                                       
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC1S3                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S4 (N, M, A, Y, H, SA, SB, SC, SD, SX, SGAMM,      &
      SALPH, SG, IERR)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 4:                                                       *      
!                                                                *      
!  All vectors, that are used for solving the systems of         *      
!  equations start with the letter S. The notation is taken      *      
!  directly from the formulas of chapter 12.1.                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRIDIG                                  *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), Y (0:M) 
      DIMENSION H (0:M - 1) 
      DIMENSION SX (M - 1), SA (M - 1), SB (M - 1), SC (M - 1), SG (M - &
      1), SGAMM (M - 1)                                                 
      DIMENSION SALPH (M - 1), SD (M - 1) 
!                                                                       
!*  determine H(J), test strict monotonicity of Y(J)                    
!                                                                       
      DO 100 J = 0, M - 1 
         H (J) = Y (J + 1) - Y (J) 
         IF (H (J) .LE.0) GOTO 900 
  100 END DO 
!                                                                       
!*  determine the columns of the tridiagonal system                     
!                                                                       
      DO 101 J = 1, M - 1 
         SB (J) = 1.0D0 / H (J - 1) 
         SC (J) = 1.0D0 / H (J) 
         SD (J) = 2.0D0 * (SB (J) + SC (J) ) 
  101 END DO 
!                                                                       
!*  in the 1st set TRIDIG is performed completely                       
!                                                                       
      IREP = 0 
!                                                                       
!*  loop covering the N+1 system of equations                           
!                                                                       
      DO 104 I = 0, N 
!                                                                       
!*  determine the right-hand side of each system                        
!                                                                       
         DO 102 J = 1, M - 1 
            SA (J) = 3.0D0 * ( (A (I, J, 1, 0) - A (I, J - 1, 1, 0) )   &
            / H (J - 1) **2 + (A (I, J + 1, 1, 0) - A (I, J, 1, 0) )    &
            / H (J) **2)                                                
  102    END DO 
!                                                                       
!*  correction in the first and last equations                          
!                                                                       
         SA (1) = SA (1) - A (I, 0, 1, 1) / H (0) 
         SA (M - 1) = SA (M - 1) - A (I, M, 1, 1) / H (M - 1) 
!                                                                       
!*  solve the system, SX is the solution vector                         
!                                                                       
         CALL TRIDIG (M - 1, SA, SB, SC, SD, IREP, SX, IERR, SGAMM,     &
         SALPH, SG)                                                     
!                                                                       
!*  return with IERR=1, if the system is unsolvable                     
!                                                                       
         IF (IERR.NE.0) RETURN 
!                                                                       
!*  to increase speed from the 2nd system onwards (compare TRIDIG)      
!                                                                       
         IREP = 1 
!                                                                       
!*  transfer solution vector to array A                                 
!                                                                       
         DO 103 J = 1, M - 1 
            A (I, J, 1, 1) = SX (J) 
  103    END DO 
  104 END DO 
!                                                                       
!*  algorithm has run successfully                                      
!                                                                       
      RETURN 
!                                                                       
!*  return with IERR=2, if monotonicity of Y(J) is not true             
!                                                                       
  900 IERR = 2 
      RETURN 
      END SUBROUTINE BIC1S4                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S5 (I, G, N, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 5:                                                       *      
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
      DIMENSION G (0:3, 0:3), X (0:N), DAT (0:3, 0:3), IPT (0:3, 0:3),  &
      HPT (0:3)                                                         
!                                                                       
!*  unchanging coefficients of the desired matrix                       
!                                                                       
      DATA DAT / 1.0D0, 0.0D0, - 3.0D0, 2.0D0, 0.0D0, 1.0D0, - 2.0D0,   &
      1.0D0, 0.0D0, 0.0D0, 3.0D0, - 2.0D0, 0.0D0, 0.0D0, - 1.0D0, 1.0D0 &
      /                                                                 
!                                                                       
!*  negative H-power of the desired matrix                              
!                                                                       
      DATA IPT / 0, 0, 2, 3, 0, 0, 1, 2, 0, 0, 2, 3, 0, 0, 1, 2 / 
!                                                                       
!*  negative H-powers                                                   
!                                                                       
      DATA HPT (0) / 1.0D0 / 
      H = X (I + 1) - X (I) 
      DO 100 K = 1, 3 
         HPT (K) = HPT (K - 1) / H 
  100 END DO 
!                                                                       
!*  determine the matrix                                                
!                                                                       
      DO 102 K = 0, 3 
         DO 101 L = 0, 3 
            IPOT = IPT (K, L) 
            G (K, L) = DAT (K, L) * HPT (IPOT) 
  101    END DO 
  102 END DO 
      RETURN 
      END SUBROUTINE BIC1S5                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S6 (J, G, M, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 6:                                                       *      
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
      DIMENSION G (0:3, 0:3), Y (0:M), DAT (0:3, 0:3), IPT (0:3, 0:3),  &
      HPT (0:3)                                                         
!                                                                       
!*  unchanging coefficients of the desired matrix                       
!                                                                       
      DATA DAT / 1.0D0, 0.0D0, - 3.0D0, 2.0D0, 0.0D0, 1.0D0, - 2.0D0,   &
      1.0D0, 0.0D0, 0.0D0, 3.0D0, - 2.0D0, 0.0D0, 0.0D0, - 1.0D0, 1.0D0 &
      /                                                                 
!                                                                       
!*  negative H-powers of the desired matrix                             
!                                                                       
      DATA IPT / 0, 0, 2, 3, 0, 0, 1, 2, 0, 0, 2, 3, 0, 0, 1, 2 / 
!                                                                       
!*  negative H-powers                                                   
!                                                                       
      DATA HPT (0) / 1.0D0 / 
      H = Y (J + 1) - Y (J) 
      DO 100 K = 1, 3 
         HPT (K) = HPT (K - 1) / H 
  100 END DO 
!                                                                       
!*  determine the transposed matrix                                     
!                                                                       
      DO 102 K = 0, 3 
         DO 101 L = 0, 3 
            IPOT = IPT (K, L) 
            G (L, K) = DAT (K, L) * HPT (IPOT) 
  101    END DO 
  102 END DO 
      RETURN 
      END SUBROUTINE BIC1S6                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S7 (N, M, A, I, J, W) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 7:                                                       *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), W (0:3, 0:3) 
      DO 102 L = 0, 1 
         DO 101 K = 0, 1 
            W (K, L) = A (I, J, K, L) 
            W (K + 2, L) = A (I + 1, J, K, L) 
            W (K, L + 2) = A (I, J + 1, K, L) 
            W (K + 2, L + 2) = A (I + 1, J + 1, K, L) 
  101    END DO 
  102 END DO 
      RETURN 
      END SUBROUTINE BIC1S7                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S8 (N, M, A, X, Y, I, J, AA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 8:                                                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BIC1S5, BIC1S6, BIC1S7                  *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), Y (0:M) 
      DIMENSION AA (0:3, 0:3), GX (0:3, 0:3), GYT (0:3, 0:3), W (0:3, 0:&
      3)                                                                
      DIMENSION WGYT (0:3, 0:3) 
!                                                                       
!*  determine G(I,J)(X) = GX                                            
!                                                                       
      CALL BIC1S5 (I, GX, N, X) 
!                                                                       
!*  determine G(I,J)(Y) transposed = GY                                 
!                                                                       
      CALL BIC1S6 (J, GYT, M, Y) 
!                                                                       
!*  determine W(I,J)                                                    
!                                                                       
      CALL BIC1S7 (N, M, A, I, J, W) 
!                                                                       
!*  W(I,J)*G(J)(Y) transposed = WGYT                                    
!                                                                       
      DO 103 K = 0, 3 
         DO 102 L = 0, 3 
            WGYT (K, L) = 0.0D0 
            DO 101 KL = 0, 3 
               WGYT (K, L) = WGYT (K, L) + W (K, KL) * GYT (KL, L) 
  101       END DO 
  102    END DO 
  103 END DO 
!                                                                       
!*  determine A(I,J) = AA                                               
!                                                                       
      DO 106 K = 0, 3 
         DO 105 L = 0, 3 
            AA (K, L) = 0.0D0 
            DO 104 KL = 0, 3 
               AA (K, L) = AA (K, L) + GX (K, KL) * WGYT (KL, L) 
  104       END DO 
  105    END DO 
  106 END DO 
      RETURN 
      END SUBROUTINE BIC1S8                         
!                                                                       
!                                                                       
      SUBROUTINE BIC1S9 (N, M, A, I, J, AA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Step 9:                                                       *      
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
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM) 
      DIMENSION AA (0:3, 0:3) 
!                                                                       
      DO 102 K = 0, 3 
         DO 101 L = 0, 3 
            A (I, J, K, L) = AA (K, L) 
  101    END DO 
  102 END DO 
      RETURN 
      END SUBROUTINE BIC1S9                         

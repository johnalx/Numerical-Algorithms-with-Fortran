      SUBROUTINE PR2TRA (NX, X, Y, F, M, MARK, C, A, IWORK, WK, XQUER,  &
      YQUER, R)                                                         
!                                                                       
!*****************************************************************      
!                                                                *      
!  PR2TRA determines a 2-dimensional surface spline for any set  *      
!  of triples (X(I),Y(I),F(X(I),Y(I))), I=1, ..., NX. The pairs  *      
!  (X(I),Y(I)) must be distinct, i.e., for each value (X,Y) there*      
!  must be exactly one value F=F(X,Y), or F must represent a     *      
!  functional relation.                                          *      
!  The nodes (X(I),Y(I)) do not have to lie on a rectangular     *      
!  grid. In fact they can be given in any order. The nodes       *      
!  (X(I),Y(I)) are transformed onto the unit circle.             *      
!  The degree of smoothness required, i.e., the order of differ- *      
!  entiability should not be chosen too high, since the condition*      
!  of the linear system of equations that will have to be solved *      
!  deteriorates with increasing derivative order. Tests have     *      
!  shown that derivative orders between 3 and 5 are advisable.   *      
!  Only in some rare cases have higher derivative orders resulted*      
!  in noticable improvements. The condition of the system of     *      
!  equations also worsens if the number of nodes increases or if *      
!  their distance decreases.                                     *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   :  number of nodes                                        *      
! X,Y  :  vectors X(1:NX), Y(1:NX); the nodes at which the       *      
!         functional values are known                            *      
! F    :  vector F(1:NX); containing the functional values at    *      
!         the nodes (X(I), Y(I))                                 *      
! M    :  derivative order with which the spline coefficients    *      
!         determined                                             *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! C    :  vector C(1:(NX + M*(M+1)/2)); the coefficients of the  *      
!         spline function                                        *      
! MARK :  indicates whether the system of equations is solvable  *      
!         MARK = 1:  everything o.k.                             *      
!         MARK = 0:  system matrix is numerically singular       *      
! X,Y  :  vectors X(1:NX), Y(1:NX); the nodes transformed to the *      
!         the unit circle                                        *      
! XQUER:  mean of the X(I) values                                *      
! YQUER:  mean of the Y(I) values                                *      
! R    :  maximal distance of a node (X(I),Y(I)) from their      *      
!         center of gravity (XQUER,YQUER)                        *      
!                                                                *      
!                                                                *      
! AUXILIARY PARAMETERS:                                          *      
! =====================                                          *      
! A    : Vector A(1:(NX + M*(M+1)/2)*(3 + NX + M*(M+1)/2)/2)     *      
! IWORK: INTEGER vector IWORK(1:(NX + M*(M+1)/2))                *      
! WK   : vector WK(1:(NX + M*(M+1)/2)*((NX + M*(M+1)/2)+1)/2)    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: ALPHA2, GAMMA2, CEPSPM, SESSPM, TRCIRC   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX), F (NX), C (NX + M * (M + 1) / 2),       &
      A ( (NX + M * (M + 1) / 2) * (3 + NX + M * (M + 1) / 2) / 2),     &
      WK ( (NX + M * (M + 1) / 2) * ( (NX + M * (M + 1) / 2) + 1)       &
      / 2)                                                              
      INTEGER IWORK (NX + M * (M + 1) / 2) 
!..                                                                     
!..   size of the system matrix                                         
!..                                                                     
      NM = NX + M * (M + 1) / 2 
!..                                                                     
!..   indicator for the polynomial components of the matrix             
!..                                                                     
      NXX = 1 + NX * (NX + 1) / 2 
!..                                                                     
!..   initializing the error parameter                                  
!..                                                                     
      MARK = 1 
!..                                                                     
!..   transformation of the nodes X(I), Y(I)                            
!..   to the unit circle                                                
!..                                                                     
      CALL TRCIRC (X, Y, NX, XQUER, YQUER, R) 
!..                                                                     
!..   formimg the matrix:                                               
!..   polynomial components P: top right, in condensed form             
!..                                                                     
      CALL ALPHA2 (NX, X, Y, M, A (NXX), IWORK (1), IWORK (1 + (M + 1)  &
      * M / 2) )                                                        
!..                                                                     
!..   core components G of the matrix:                                  
!..   top left, upper triangle, in condensed form                       
!..                                                                     
      CALL GAMMA2 (NX, X, Y, M, A) 
!..                                                                     
!..   initializing the right-hand side                                  
!..                                                                     
      DO 20 I = 1, NX 
         C (I) = F (I) 
   20 END DO 
      DO 30 I = NX + 1, NM 
         C (I) = 0.0D0 
   30 END DO 
!..                                                                     
!..   decomposing the system matrix                                     
!..                                                                     
      CALL CEPSPM (A, NM, C, IWORK, RCOND, A (NM * (NM + 1) / 2 + 1),   &
      WK)                                                               
!..                                                                     
!..   if the system matrix numerically is singular: stop                
!..                                                                     
      IF (1.0D0.EQ.1.0D0 + RCOND) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!..                                                                     
!..   solve the system of equations                                     
!..                                                                     
      CALL SESSPM (WK, NM, IWORK, C) 
      RETURN 
      END SUBROUTINE PR2TRA                         
!                                                                       
!                                                                       
      SUBROUTINE APPRT2 (X0, Y0, NX, M, X, Y, C, AP, XQUER, YQUER, R) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Evaluation function used for interpolation of the surface      *      
! spline.                                                        *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X0,Y0:  location where the surface spline is to be evaluated   *      
! NX   :  number of nodes                                        *      
! M    :  derivative order for which the spline coefficients are *      
!         determined                                             *      
! X,Y  :  vectors X(1:NX), Y(1:NX); the transformed nodes on the *      
!         unit circle for which the coefficients were determined *      
! C    :  vector C(1:(NX + M*(M+1)/2)); coefficient vector       *      
! XQUER:  arithmetic mean of the X(I)                            *      
! YQUER:  arithmetic mean of the Y(I)                            *      
! R    :  largest distance of a node (X(I),Y(I)) from the center *      
!         (XQUER,YQUER)                                          *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! AP   :  approximation value of the spline at (X0,Y0)           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: none                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Richard Reuter, 1983                                *      
! editor   : Hartmut Turowski                                    *      
! date     : 06.10.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX), C (NX + M * (M + 1) / 2) 
!..                                                                     
!..   transformation of the point (X0,Y0), where the evaluation is      
!..   to be performed, to the unit circle                               
!..                                                                     
      X0 = R * (X0 - XQUER) 
      Y0 = R * (Y0 - YQUER) 
!..                                                                     
!..   for different M various cases are considered separately:          
!..   1. M = 1, 2, 3 ; special coding, very fast                        
!..   2. M > 3 ; each monomial is represented in the form (X**IX)*(Y**IY
!..              A functional evaluation is slow and prone to rounding  
!..              errors                                                 
!..                                                                     
!..   the starting polynomial always is equal to 1                      
!..                                                                     
      AP = C (NX + 1) 
      IF (M.EQ.1) GOTO 20 
      IF (M.EQ.2) THEN 
         AP = AP + C (NX + 2) * X0 + C (NX + 3) * Y0 
      ELSEIF (M.EQ.3) THEN 
         AP = AP + (C (NX + 2) + C (NX + 4) * X0 + C (NX + 5) * Y0)     &
         * X0 + (C (NX + 3) + C (NX + 6) * Y0) * Y0                     
      ELSE 
         IX = 0 
         IY = 0 
         DO 10 I = 2, M * (M + 1) / 2 
            IF (IX.EQ.0) THEN 
               IX = IY + 1 
               IY = 0 
               AP = AP + C (NX + I) * (X0**IX) 
            ELSE 
               IX = IX - 1 
               IY = IY + 1 
               IF (IX.EQ.0) THEN 
                  AP = AP + C (NX + I) * (Y0**IY) 
               ELSE 
                  AP = AP + C (NX + I) * (X0**IX) * (Y0**IY) 
               ENDIF 
            ENDIF 
   10    END DO 
      ENDIF 
   20 CONTINUE 
!..                                                                     
!..   component of core function E                                      
!..                                                                     
!..   The function E2(X,Y,M) could be called at this point.             
!..   However, this would slow down the evaluation con-                 
!..   siderably. Thus, a direct code is performed                       
!..                                                                     
      DO 30 I = 1, NX 
         R2 = (X (I) - X0) **2 + (Y (I) - Y0) **2 
         IF (R2.EQ.0.0D0) R2 = 1.0D0 
         AP = AP + C (I) * DLOG (R2) * R2** (M - 1) 
   30 END DO 
      RETURN 
      END SUBROUTINE APPRT2                         
!                                                                       
!                                                                       
      SUBROUTINE TRCIRC (X, Y, NX, XQUER, YQUER, R) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Transformation of the nodes (X(I),Y(I)) to the unit circle.    *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX   : number of nodes                                         *      
! X,Y  : vectors X(1:NX), X(1:NX); the nodes where the spline    *      
!        function is known                                       *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X,Y  :  vectors X(1:NX), Y(1:NX); nodes that have been trans-  *      
!         formed to the unit circle                              *      
! XQUER:  arithmetic mean of the X(I)                            *      
! YQUER:  arithmetic mean of the Y(I)                            *      
! R    :  largest distance of a node (X(I),Y(I)) from the center *      
!         (XQUER,YQUER)                                          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! subroutines required: none                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! author   : Hartmut Turowski                                    *      
! date     : 07.23.1988                                          *      
! source   : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DIMENSION X (NX), Y (NX) 
!..                                                                     
!..   determine the arithmetic means for                                
!..   the X(I) and Y(I)                                                 
!..                                                                     
      XQUER = 0.0D0 
      YQUER = 0.0D0 
      DO 10 I = 1, NX 
         XQUER = XQUER + X (I) 
         YQUER = YQUER + Y (I) 
   10 END DO 
      XQUER = XQUER / FLOAT (NX) 
      YQUER = YQUER / FLOAT (NX) 
!..                                                                     
!..   determine the maximal distance of the nodes                       
!..   (X(I),Y(I)) from their center of gravity (XQUER,YQUER)            
!..                                                                     
      R = 0.0D0 
      DO 20 I = 1, NX 
         R = DMAX1 (DSQRT ( (X (I) - XQUER) **2 + (Y (I) - YQUER) **2), &
         R)                                                             
   20 END DO 
      R = 1.0D0 / R 
!..                                                                     
!..   transformation to the unit circle                                 
!..                                                                     
      DO 30 I = 1, NX 
         X (I) = R * (X (I) - XQUER) 
         Y (I) = R * (Y (I) - YQUER) 
   30 END DO 
      RETURN 
      END SUBROUTINE TRCIRC                         

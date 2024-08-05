      SUBROUTINE PR3TRA (NX, X, Y, Z, F, M, MARK, C, A, IWORK, WK,      &
      XMEAN, YMEAN, ZMEAN, R)                                           
!                                                                       
!*****************************************************************      
!                                                                *      
! PR3TRA computes three-dimensional surface splines for arbitrary*      
! given points (X(I),Y(I),Z(I),F(X(I),Y(I)), I=1, ..., NX.       *      
! The nodes (X(I),Y(I),Z(I)) must be distinct and F must be a    *      
! function, i.e., for each (X,Y,Z) in the node set there must    *      
! correspond a unique  F=F(X,Y,Z). The nodes need not be ordered.*      
! The desired smoothness of the spline, i. e., its derivative    *      
! order should be stipulated as rather low since the condition   *      
! number of the system of equations that has to be solved worsens*      
! with increasing derivative order.                              *      
! Tests indicate that derivative orders between  3 and 5  can be *      
! recommended. Higher orders showed improvement only in rare     *      
! cases. For an increasing number of nodes, i. e., a decreasing  *      
! distance between the nodes the condition number of the linear  *      
! system also tends to worsen.                                   *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX     :  Number of nodes                                      *      
! X,Y,Z  :  NX-vectors ..(1:NX); the coordinates of the nodes    *      
! F      :  NX-vector F(1:NX); the functional values at the nodes*      
! M      :  derivative order used to determine the coefficients  *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
!                                                                *      
! C      :  vector C(1:(NX + M*(M+1)*(M+2)/6)); the coefficients *      
!           of the spline                                        *      
!                                                                *      
! MARK   :  indicates whether the linear system could be solved: *      
!           MARK = 1:  all is ok                                 *      
!           MARK = 0:  system matrix is numerically singular     *      
! XMEAN  : )                                                     *      
! YMEAN  : )  refer to  SUBROUTINE TRBALL                        *      
! ZMEAN  : )                                                     *      
! R      : )                                                     *      
!                                                                *      
! AUXILIARY PARAMETERS:                                          *      
! =====================                                          *      
! A      : vector A(1:((NX + M*(M+1)*(M+2)/6) *                  *      
!                         * (NX + M*(M+1)*(M+2)/6 + 3))/2)       *      
! IWORK  : integer vector IWORK(1:(NX + M*(M+1)*(M+2)/6))        *      
! WK     : vector WK(1:((NX + M*(M+1)*(M+2)/6)*                  *      
!                         * (NX + M*(M+1)*(M+2)/6 + 1)/2))       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines: ALPHA3, GAMMA3, NEXT3, E3, CEPSPM,       *      
!                       ZSPMMK, PCOSOL, PCOLTG, SESSPM, SCAPRO,  *      
!                       VECMWC, ABSSUM, INDMAX, VECADD, VECXCH,  *      
!                       TRBALL                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (NX), Y (NX), Z (NX), F (NX) 
      DIMENSION C (NX + M * (M + 1) * (M + 2) / 6) 
      DIMENSION A ( (NX + M * (M + 1) * (M + 2) / 6) * (NX + M *        &
      (M + 1) * (M + 2) / 6 + 3) / 2)                                   
      DIMENSION WK ( (NX + M * (M + 1) * (M + 2) / 6) * (NX + M *       &
      (M + 1) * (M + 2) / 6 + 1) / 2)                                   
      DIMENSION IWORK (NX + M * (M + 1) * (M + 2) / 6) 
!..                                                                     
!..   size of the matrix                                                
!..                                                                     
      M3 = M * (M + 1) * (M + 2) / 6 
      NM = NX + M3 
!..                                                                     
!..   pointer for the polynomial part of the matrix                     
!..                                                                     
      NXX = 1 + (NX * (NX + 1) ) / 2 
!..                                                                     
!..   Initialize error parameter                                        
!..                                                                     
      MARK = 1 
!..                                                                     
!..   transform all nodes X(I), Y(I), Z(I)                              
!..   onto the unit ball                                                
!..                                                                     
      CALL TRBALL (X, Y, Z, NX, XMEAN, YMEAN, ZMEAN, R) 
!..                                                                     
!..   Form the system matrix:                                           
!..   its polynomial part P appears in the upper right corner           
!..   in condensed form                                                 
!..                                                                     
      CALL ALPHA3 (NX, X, Y, Z, M - 1, A (NXX), IWORK (1), IWORK (1 +   &
      M3), IWORK (1 + 2 * M3) )                                         
!..                                                                     
!..   G part of the system matrix:                                      
!..   in the upper left corner in condensed form                        
!..                                                                     
      CALL GAMMA3 (NX, X, Y, Z, M, A) 
!..                                                                     
!..   set up the right hand side                                        
!..                                                                     
      DO 10 I = 1, NX 
         C (I) = F (I) 
   10 END DO 
      DO 20 I = NX + 1, NM 
         C (I) = 0.0D0 
   20 END DO 
!..                                                                     
!..   factor the system matrix                                          
!..                                                                     
      CALL CEPSPM (A, NM, C, IWORK, RCOND, A ( (NM * (NM + 1) ) / 2 + 1)&
      , WK)                                                             
!..                                                                     
!..   Stop if the system matrix is numerically singular                 
!..                                                                     
      IF (1.0D0.EQ.1.0D0 + RCOND) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!..                                                                     
!..   solve the system                                                  
!..                                                                     
      CALL SESSPM (WK, NM, IWORK, C) 
      RETURN 
      END SUBROUTINE PR3TRA                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION APPTR3 (X0, Y0, Z0, NX, M, X, Y, Z, C,   &
      XMEAN, YMEAN, ZMEAN, R)                                           
!                                                                       
!*****************************************************************      
!                                                                *      
! Evaluation function for the interpolation                      *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X0,Y0,Z0 :  location where the function is to be evaluated     *      
! NX       :  number of nodes                                    *      
! M        :  derivative order                                   *      
! X,Y,Z    :  NX-vectors ..(1:NX); the coordinates of the nodes  *      
! C        :  vector C(1:(NX + M*(M+1)*(M+2)/6)) of coefficients *      
! XMEAN  : )                                                     *      
! YMEAN  : )  refer to  SUBROUTINE TRBALL                        *      
! ZMEAN  : )                                                     *      
! R      : )                                                     *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! APPTR3   :  approximate value at (X0,Y0,Z0)                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines:  none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Authors     : Richard Reuter (1983), Hartmut Turowski          *      
! Date        : 12.10.1989                                       *      
! Source      : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (NX), Y (NX), Z (NX), C (NX + (M * (M + 1) * (M + 2)  &
      / 6) )                                                            
!..                                                                     
!..   Transform the point (X0,Y0,Z0), where we want to                  
!..   evaluate to the unit ball                                         
!..                                                                     
      X0 = R * (X0 - XMEAN) 
      Y0 = R * (Y0 - YMEAN) 
      Z0 = R * (Z0 - ZMEAN) 
                                                                        
!..                                                                     
!..   Depending on m there are several cases:                           
!..   1. M = 1, 2, 3 ; separately coded, fast                           
!..   2. M > 3       ; each monomial is represented in the form         
!..                    (X**IX)*(Y**IY)*(Z**IZ).                         
!..                    The evaluation is slower and prone to            
!..                    rounding errors.                                 
!..                                                                     
!..   the first polynomial is  1                                        
!..                                                                     
      AP = C (NX + 1) 
      IF (M.EQ.1) GOTO 40 
      IF (M.EQ.2) THEN 
!..                                                                     
!..   remaining monomials of degree 1                                   
!..                                                                     
         AP = AP + C (NX + 2) * X0 + C (NX + 3) * Y0 + C (NX + 4)       &
         * Z0                                                           
      ELSEIF (M.EQ.3) THEN 
!..                                                                     
!..   remaining monomials of degree 2                                   
!..                                                                     
         AP = AP + (C (NX + 2) + C (NX + 5) * X0 + C (NX + 6) * Y0)     &
         * X0 + (C (NX + 3) + C (NX + 8) * Y0 + C (NX + 9) * Z0)        &
         * Y0 + (C (NX + 4) + C (NX + 7) * X0 + C (NX + 10) * Z0)       &
         * Z0                                                           
      ELSE 
!..                                                                     
!..   remaining monomials of degree <= M-1                              
!..                                                                     
         L = 1 
         DO 30 I = 1, M - 1 
            DO 20 IX = I, 0, - 1 
               DO 10 IY = I - IX, 0, - 1 
                  IZ = I - IX - IY 
                  L = L + 1 
                  IF (IX.NE.0.AND.IY.NE.0.AND.IZ.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) * (Y0**IY) *       &
                     (Z0**IZ)                                           
                  ELSEIF (IX.NE.0.AND.IY.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) * (Y0**IY) 
                  ELSEIF (IX.NE.0.AND.IZ.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) * (Z0**IZ) 
                  ELSEIF (IY.NE.0.AND.IZ.NE.0) THEN 
                     AP = AP + C (NX + L) * (Y0**IY) * (Z0**IZ) 
                  ELSEIF (IX.NE.0) THEN 
                     AP = AP + C (NX + L) * (X0**IX) 
                  ELSEIF (IY.NE.0) THEN 
                     AP = AP + C (NX + L) * (Y0**IY) 
                  ELSE 
                     AP = AP + C (NX + L) * (Z0**IZ) 
                  ENDIF 
   10          END DO 
   20       END DO 
   30    END DO 
      ENDIF 
   40 CONTINUE 
!..                                                                     
!..   Kernelfunction                                                    
!..                                                                     
!..   one could call the FUNCTION E3(X,Y,Z) now, but this would         
!..   slow down the actual evaluation, hence we code directly here.     
!..                                                                     
      DO 50 I = 1, NX 
         R = (X (I) - X0) **2 + (Y (I) - Y0) **2 + (Z (I) - Z0) **2 
         AP = AP + C (I) * (DSQRT (R) ** (2 * M - 3) ) 
   50 END DO 
      APPTR3 = AP 
      RETURN 
      END FUNCTION APPTR3                           
!                                                                       
!                                                                       
      SUBROUTINE TRBALL (X, Y, Z, NX, XMEAN, YMEAN, ZMEAN, R) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Transform all nodes (X(I),Y(I),Z(I)) to the unit ball.         *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! NX    : number of nodes                                        *      
! X,Y,Z : NX-vectors ..(1:NX); the coordinates of the nodes      *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X,Y,Z : NX-vectors ..(1:NX); the nodes transformed to the unit *      
!         ball                                                   *      
! XMEAN : arithmetic mean of the X(I)                            *      
! YMEAN : arithmetic mean of the Y(I)                            *      
! ZMEAN : arithmetic mean of the Z(I)                            *      
! R     : largest euclidean distance of a node (X(I),Y(I),Z(I)   *      
!         from the center of gravity (XMEAN,YMEAN,ZMEAN)         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
! Required subroutines: none                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
! Author     : Hartmut Turowski                                  *      
! Date       : 01.01.1990                                        *      
! Source     : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DOUBLEPRECISION X (NX), Y (NX), Z (NX), XMEAN, YMEAN, ZMEAN, R,   &
      RHILF                                                             
!..                                                                     
!..   compute the arithmetic means of the                               
!..   X(I), Y(I)  and  Z(I)                                             
!..                                                                     
      XMEAN = 0.0D0 
      YMEAN = 0.0D0 
      ZMEAN = 0.0D0 
      DO 10 I = 1, NX 
         XMEAN = XMEAN + X (I) 
         YMEAN = YMEAN + Y (I) 
         ZMEAN = ZMEAN + Z (I) 
   10 END DO 
      XMEAN = XMEAN / DBLE (NX) 
      YMEAN = YMEAN / DBLE (NX) 
      ZMEAN = ZMEAN / DBLE (NX) 
!..                                                                     
!..   compute the maximal distance of a node (X(I),Y(I),Z(I))           
!..   from (XMEAN,YMEAN,ZMEAN)                                          
!..                                                                     
      R = 0.0D0 
      DO 20 I = 1, NX 
         RHILF = DSQRT ( (X (I) - XMEAN) **2 + (Y (I) - YMEAN) **2 +    &
         (Z (I) - ZMEAN) **2)                                           
         R = DMAX1 (RHILF, R) 
   20 END DO 
      R = 1.0D0 / R 
!..                                                                     
!..   Transform to the unit ball                                        
!..                                                                     
      DO 30 I = 1, NX 
         X (I) = R * (X (I) - XMEAN) 
         Y (I) = R * (Y (I) - YMEAN) 
         Z (I) = R * (Z (I) - ZMEAN) 
   30 END DO 
      RETURN 
      END SUBROUTINE TRBALL                         

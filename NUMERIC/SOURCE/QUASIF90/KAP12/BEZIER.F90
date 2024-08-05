![            {B\'ezier Spline Surfaces}*)                              
      SUBROUTINE BEZIER (B, D, WORK, ITYP, M, N, EPS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This subroutine performs                                   *      
!     for ITYP = 0   the bicubic BEZIER-method. Here we determine*      
!                    interpolation points from the input data for*      
!                    a spline surface that is constructed using  *      
!                    the bicubic BEZIER-method;                  *      
!     for ITYP = 1   the modified bicubic BEZIER-method. In it   *      
!                    the input interpolation points are regarded *      
!                    as weight points at first, for which pseudo-*      
!                    interpolation points are determined. These  *      
!                    are shifted until they coincide with the    *      
!                    true interpolation points, except for an    *      
!                    accuracy EPS.                               *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     ITYP          INTEGER   which defines the method used:     *      
!                               ITYP=0 : BEZIER-method           *      
!                               ITYP=1 : modified BEZIER-method  *      
!     M             INTEGER   number of patches in 1st direction *      
!     N             INTEGER   number of patches in 2nd direction *      
!                                                                *      
!     ---   for ITYP = 0  ---                                    *      
!     D(3,0:M,0:N)     Double Precision coordinates of the weight*      
!                      points                                    *      
!                                                                *      
!     ---  for ITYP = 1  ---                                     *      
!     EPS              Double Precision accuracy bound for the   *      
!                             interpolation                      *      
!     D(3,0:M,0:N)     Double Precision coordinates of the inter-*      
!                             polation points                    *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B(3,0:3M,0:3N)   Double Precision coordinates of the BEZIER*      
!                      points  B(3,J,I) for J=0, 1, ..., 3M and  *      
!                      I=0, 1, ..., 3N                           *      
!                                                                *      
!                                                                *      
!     AUXILIARY PARAMETERS:                                      *      
!     =====================                                      *      
!     WORK(3,0:M,0:N)  Double Precision storage for coordinates  *      
!                      of intermediate BEZIER points.            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: INTPOL, BEZPNT, BEZBRD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Michael Radermacher                                *      
!  editor   : Hartmut Turowski                                   *      
!  date     : 05.24.1990                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION B (3, 0:3 * M, 0:3 * N), D (3, 0:M, 0:N), WORK (3,&
      0:M, 0:N)                                                         
      DOUBLEPRECISION DIFF (3) 
!                                                                       
!*****************************************************************      
!     for the bicubic BEZIER-method:                             *      
!     call SUBROUTINE BEZPNT to determine the BEZIER points      *      
!*****************************************************************      
!                                                                       
      IF (ITYP.EQ.0) THEN 
         CALL BEZPNT (B, D, M, N) 
         RETURN 
      ENDIF 
!                                                                       
!*****************************************************************      
!     for the modified bicubic BEZIER-method:                    *      
!     save the contents of D in work-space WORK                  *      
!*****************************************************************      
!                                                                       
      DO 20 I = 0, N 
         DO 20 J = 0, M 
            DO 20 L = 1, 3 
   20 WORK (L, J, I) = D (L, J, I) 
!                                                                       
!*****************************************************************      
!     call SUBROUTINE BEZPNT for first determination of          *      
!     the BEZIER points using the bicubic BEZIER-method          *      
!*****************************************************************      
!                                                                       
      CALL BEZPNT (B, D, M, N) 
!                                                                       
!*****************************************************************      
!     transformation of the coordinates of the just determined   *      
!     BEZIER points to the original interpolation points         *      
!*****************************************************************      
!                                                                       
   40 DO 60 J = 0, M 
         DO 60 I = 0, N 
            DO 50 L = 1, 3 
   50       DIFF (L) = WORK (L, J, I) - B (L, 3 * J, 3 * I) 
   60 CALL INTPOL (DIFF, 3 * J, 3 * I, B, M, N) 
!                                                                       
!*****************************************************************      
!     test whether the modified interpolation points approximate *      
!     the original interpolation points within EPS               *      
!*****************************************************************      
!                                                                       
      DO 70 J = 0, M 
         DO 70 I = 0, N 
            DO 70 L = 1, 3 
   70 IF (ABS (B (L, 3 * J, 3 * I) - WORK (L, J, I) ) .GT.EPS) GOTO 40 
      RETURN 
      END SUBROUTINE BEZIER                         
!                                                                       
!                                                                       
      SUBROUTINE INTPOL (DIFF, J, I, B, M, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This Subroutine performs the changes of the interpolation  *      
!     points stored in B for the spline surface that has been    *      
!     approximately determined by the bicubic BEZIER-method      *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     DIFF(3)         Double Precision coordinates of the        *      
!                     difference vector, according to which the  *      
!                     BEZIER surface is to be modified           *      
!     J, I            INTEGER labels for the patch, in whose     *      
!                     neighborhood the BEZIER surface is modified*      
!     M               INTEGER: number of patches in 1st direction*      
!     N               INTEGER: number of patches in 2nd direction*      
!     B(3,0:3M,0:3N)  Double Precision coordinates of the BEZIER *      
!                     points                                     *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B(3,0:3M,0:3N)  Double Precision coordinates of the BEZIER *      
!                     points                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Michael Radermacher                                *      
!  editor   : Hartmut Turowski                                   *      
!  date     : 05.18.1990                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION B (3, 0:3 * M, 0:3 * N), DIFF (3) 
!                                                                       
!*****************************************************************      
!     in case the BEZIER point B(L,3*J,3*I) that is to be        *      
!     modified lies on then boundary of the surface              *      
!*****************************************************************      
!                                                                       
      IF ( (I.EQ.0) .OR. (I.EQ.3 * N) ) THEN 
!                                                                       
!*****************************************************************      
!     the point lies on the boundary with I=0 or I=3*N           *      
!*****************************************************************      
!                                                                       
         DO 40 L = 1, 3 
            IF ( (J.EQ.0) .OR. (J.EQ.3 * M) ) THEN 
!                                                                       
!*****************************************************************      
!     corner point: modify by 16*D                               *      
!*****************************************************************      
!                                                                       
               B (L, J, I) = B (L, J, I) + DIFF (L) 
            ELSE 
!                                                                       
!*****************************************************************      
!     other boundary point:                                      *      
!     modify the points with J-3 and J+3 by 6*D                  *      
!*****************************************************************      
!                                                                       
               DO 10 K1 = - 3, 3, 6 
   10          B (L, J + K1, I) = B (L, J + K1, I) + 3.0D0 * DIFF (L)   &
               / 8.0D0                                                  
!                                                                       
!*****************************************************************      
!     modify the points with J-2 and J+2 by 12*D                 *      
!*****************************************************************      
!                                                                       
               DO 20 K1 = - 2, 2, 4 
   20          B (L, J + K1, I) = B (L, J + K1, I) + .75D0 * DIFF (L) 
!                                                                       
!*****************************************************************      
!     modify the 3 interior points by 24*D                       *      
!*****************************************************************      
!                                                                       
               DO 30 K1 = - 1, 1 
   30          B (L, J + K1, I) = B (L, J + K1, I) + 1.5D0 * DIFF (L) 
            ENDIF 
   40    END DO 
         RETURN 
      ENDIF 
      IF ( (J.EQ.0) .OR. (J.EQ.3 * M) ) THEN 
!                                                                       
!*****************************************************************      
!     the point lies on the boundary with J=0 or J=3*M           *      
!     modify analogously as for the I-margins:                   *      
!*****************************************************************      
!                                                                       
         DO 80 L = 1, 3 
            IF ( (I.EQ.0) .OR. (I.EQ.3 * N) ) THEN 
!                                                                       
!*****************************************************************      
!     corner point                                               *      
!*****************************************************************      
!                                                                       
               B (L, J, I) = B (L, J, I) + DIFF (L) 
            ELSE 
!                                                                       
!*****************************************************************      
!     other boundary point                                       *      
!*****************************************************************      
!                                                                       
               DO 50 K2 = - 3, 3, 6 
   50          B (L, J, I + K2) = B (L, J, I + K2) + 3.0D0 * DIFF (L)   &
               / 8.0D0                                                  
               DO 60 K2 = - 2, 2, 4 
   60          B (L, J, I + K2) = B (L, J, I + K2) + .75D0 * DIFF (L) 
               DO 70 K2 = - 1, 1 
   70          B (L, J, I + K2) = B (L, J, I + K2) + 1.5D0 * DIFF (L) 
            ENDIF 
   80    END DO 
         RETURN 
      ENDIF 
!                                                                       
!*****************************************************************      
!     loop over the X-, Y- and Z-coordinates                     *      
!*****************************************************************      
!                                                                       
      DO 180 L = 1, 3 
!                                                                       
!*****************************************************************      
!        modification of the BEZIER points B(L,3*J-1,3*I-1),     *      
!        B(L,3*J-1,3*I+1), B(L,3*J,3*I-1), B(L,3*J,3*I),         *      
!        B(L,3*J+1,3*I-1), B(L,3*J+1,3*I), B(L,3*J+1,3*I+1),     *      
!        B(L,3*J-1,3*I) and B(L,3*J,3*I+1) by 16*D               *      
!*****************************************************************      
!                                                                       
         DO 90 K1 = - 1, 1 
            DO 90 K2 = - 1, 1 
   90    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L) 
!                                                                       
!*****************************************************************      
!      modification of the BEZIER points B(L,3*J-2,3*I-1),       *      
!      B(L,3*J-2,3*I), B(L,3*J-2,3*I+1), B(L,3*J+2,3*I-1),       *      
!      B(L,3*J+2,3*I) and B(L,3*J+2,3*I+1) by 8*D                *      
!*****************************************************************      
!                                                                       
         DO 100 K1 = - 2, 2, 4 
            DO 100 K2 = - 1, 1 
  100    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 2.0D0                                                        
!                                                                       
!*****************************************************************      
!      modification of the BEZIER points B(L,3*J-1,3*I-2),       *      
!      B(L,3*J,3*I-2), B(L,3*J+1,3*I-2), B(L,3*J-1,3*I+2),       *      
!      B(L,3*J,3*I+2), and B(L,3*J+1,3*I+2) by 8*D               *      
!*****************************************************************      
!                                                                       
         DO 110 K1 = - 1, 1 
            DO 110 K2 = - 2, 2, 4 
  110    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 2.0D0                                                        
!                                                                       
!*****************************************************************      
!     modification of the BEZIER points B(L,3*J-3,3*I-1),        *      
!     B(L,3*J-3,3*I), B(L,3*J-3,3*I+1), B(L,3*J+3,3*I-1),        *      
!     B(L,3*J+3,3*I) and B(L,3*J+3,3*I+1) by 4*D                 *      
!*****************************************************************      
!                                                                       
         DO 120 K1 = - 3, 3, 6 
            DO 120 K2 = - 1, 1 
  120    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 4.0D0                                                        
!                                                                       
!*****************************************************************      
!   modification of the BEZIER points B(L,3*J-2,3*I-2),          *      
!   B(L,3*J-2,3*I+2), B(L,3*J+2,3*I-2) and B(L,3*J+2,3*I+2)      *      
!   by 4*D                                                       *      
!*****************************************************************      
!                                                                       
         DO 130 K1 = - 2, 2, 4 
            DO 130 K2 = - 2, 2, 4 
  130    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 4.0D0                                                        
!                                                                       
!*****************************************************************      
!   modification of the BEZIER points B(L,3*J-1,3*I-3),          *      
!   B(L,3*J-1,3*I+3), B(L,3*J,3*I-3), B(L,3*J,3*I+3),            *      
!   B(L,3*J+1,3*I-3) and B(L,3*J+1,3*I+3) by 4*D                 *      
!*****************************************************************      
!                                                                       
         DO 140 K1 = - 1, 1 
            DO 140 K2 = - 3, 3, 6 
  140    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 4.0D0                                                        
!                                                                       
!*****************************************************************      
!   modification of the BEZIER points B(L,3*J-3,3*I-2),          *      
!   B(L,3*J-3,3*I+2), B(L,3*J+3,3*I-2) and B(L,3*J+3,3*I+2)      *      
!   by 2*D                                                       *      
!*****************************************************************      
!                                                                       
         DO 150 K1 = - 3, 3, 6 
            DO 150 K2 = - 2, 2, 4 
  150    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 8.0D0                                                        
!                                                                       
!*****************************************************************      
!   modification of the BEZIER points B(L,3*J-2,3*I-3),          *      
!   B(L,3*J-2,3*I+3), B(L,3*J+2,3*I-3) and B(L,3*J+2,3*I+30      *      
!   by 2*D                                                       *      
!*****************************************************************      
!                                                                       
         DO 160 K1 = - 2, 2, 4 
            DO 160 K2 = - 3, 3, 6 
  160    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 8.0D0                                                        
!                                                                       
!*****************************************************************      
!   modification of the BEZIER points B(L,3*J-3,3*I-3),          *      
!   B(L,3*J-3,3*I+3), B(L,3*J+3,3*I-3) and B(L,3*J+3,3*I+3) by D *      
!*****************************************************************      
!                                                                       
         DO 170 K1 = - 3, 3, 6 
            DO 170 K2 = - 3, 3, 6 
  170    B (L, J + K1, I + K2) = B (L, J + K1, I + K2) + DIFF (L)       &
         / 16.0D0                                                       
  180 END DO 
      RETURN 
      END SUBROUTINE INTPOL                         
!                                                                       
!                                                                       
      SUBROUTINE BEZPNT (B, D, M, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This Subroutine computes the BEZIER points, that are       *      
!     needed for determining the surface when using the bicubic  *      
!     BEZIER-method.                                             *      
!     The coordinates of the BEZIER points are stored in B;      *      
!     the corresponding weight points are to be found in D.      *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     M             INTEGER   number of patches in 1st direction *      
!     N             INTEGER   number of patches in 2nd direction *      
!     D(3,0:M,0:N)  Double Precision coordinates of the weight   *      
!                   points                                       *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B(3,0:3M,0:3N) Double Precision coordinates of all BEZIER  *      
!                    points                                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BEZBRD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Michael Radermacher                                *      
!  editor   : Hartmut Turowski                                   *      
!  date     : 05.18.1990                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION B (3, 0:3 * M, 0:3 * N), D (3, 0:M, 0:N) 
!                                                                       
!*****************************************************************      
!     determine boundary points by the 1-dimensional             *      
!     BEZIER-method                                              *      
!*****************************************************************      
!                                                                       
      CALL BEZBRD (B, D, M, N) 
!                                                                       
!*****************************************************************      
!     loop over the X-, Y- and Z-coordinates                     *      
!*****************************************************************      
!                                                                       
      DO 100 K = 1, 3 
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J-2,3*I-2);           *      
!        J=1, ..., M, I=1, ..., N, K=1,2,3                       *      
!*****************************************************************      
!                                                                       
         DO 10 J = 1, M 
            DO 10 I = 1, N 
   10    B (K, 3 * J - 2, 3 * I - 2) = (4.0D0 * D (K, J - 1, I - 1)     &
         + 2.0D0 * D (K, J - 1, I) + 2.0D0 * D (K, J, I - 1) + D (K, J, &
         I) ) / 9.0D0                                                   
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J+2,3*I-2);           *      
!        J=0, 1, ..., M-1, I=1, ..., N, K=1,2,3                  *      
!*****************************************************************      
!                                                                       
         DO 20 J = 0, M - 1 
            DO 20 I = 1, N 
   20    B (K, 3 * J + 2, 3 * I - 2) = (4.0D0 * D (K, J + 1, I - 1)     &
         + 2.0D0 * D (K, J, I - 1) + 2.0D0 * D (K, J + 1, I) + D (K, J, &
         I) ) / 9.0D0                                                   
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J-2,3*I+2);           *      
!        J=1, ..., M, I=0, ..., N-1, K=1,2,3                     *      
!*****************************************************************      
!                                                                       
         DO 30 J = 1, M 
            DO 30 I = 0, N - 1 
   30    B (K, 3 * J - 2, 3 * I + 2) = (4.0D0 * D (K, J - 1, I + 1)     &
         + 2.0D0 * D (K, J - 1, I) + 2.0D0 * D (K, J, I + 1) + D (K, J, &
         I) ) / 9.0D0                                                   
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J+2,3*I+2);           *      
!        J=0, 1, ..., M-1, I=0, 1, ..., N-1, K=1,2,3             *      
!*****************************************************************      
!                                                                       
         DO 40 J = 0, M - 1 
            DO 40 I = 0, N - 1 
   40    B (K, 3 * J + 2, 3 * I + 2) = (4.0D0 * D (K, J + 1, I + 1)     &
         + 2.0D0 * D (K, J, I + 1) + 2.0D0 * D (K, J + 1, I) + D (K, J, &
         I) ) / 9.0D0                                                   
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J-2,3*I);             *      
!        J=1, ..., N-1, I=1, ..., M, K=1,2,3                     *      
!*****************************************************************      
!                                                                       
         DO 50 J = 1, M 
            DO 50 I = 1, N - 1 
   50    B (K, 3 * J - 2, 3 * I) = (2.0D0 * D (K, J - 1, I - 1) + 8.0D0 &
         * D (K, J - 1, I) + D (K, J, I - 1) + 2.0D0 * D (K, J - 1, I + &
         1) + 4.0D0 * D (K, J, I) + D (K, J, I + 1) ) / 18.0D0          
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J,3*I-2);             *      
!        J=1, ..., M-1, I=1, ..., N, K=1,2,3                     *      
!*****************************************************************      
!                                                                       
         DO 60 J = 1, M - 1 
            DO 60 I = 1, N 
   60    B (K, 3 * J, 3 * I - 2) = (2.0D0 * D (K, J - 1, I - 1) + 8.0D0 &
         * D (K, J, I - 1) + D (K, J - 1, I) + 2.0D0 * D (K, J + 1, I - &
         1) + 4.0D0 * D (K, J, I) + D (K, J + 1, I) ) / 18.0D0          
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J,3*I+2);             *      
!        J=1, ..., M-1, I=0, 1, ..., N-1, K=1,2,3                *      
!*****************************************************************      
!                                                                       
         DO 70 J = 1, M - 1 
            DO 70 I = 0, N - 1 
   70    B (K, 3 * J, 3 * I + 2) = (2.0D0 * D (K, J - 1, I + 1) + 8.0D0 &
         * D (K, J, I + 1) + D (K, J - 1, I) + 2.0D0 * D (K, J + 1, I + &
         1) + 4.0D0 * D (K, J, I) + D (K, J + 1, I) ) / 18.0D0          
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J+2,3*I);             *      
!        J=0, 1, ..., M-1, I=1, ..., N-1, K=1,2,3                *      
!*****************************************************************      
!                                                                       
         DO 80 J = 0, M - 1 
            DO 80 I = 1, N - 1 
   80    B (K, 3 * J + 2, 3 * I) = (2.0D0 * D (K, J + 1, I - 1) + 8.0D0 &
         * D (K, J + 1, I) + D (K, J, I - 1) + 2.0D0 * D (K, J + 1, I + &
         1) + 4.0D0 * D (K, J, I) + D (K, J, I + 1) ) / 18.0D0          
!                                                                       
!*****************************************************************      
!        determine the BEZIER points B(K,3*J,3*I);               *      
!        J=1, ..., M-1, I=1, ..., N-1, K=1,2,3                   *      
!*****************************************************************      
!                                                                       
         DO 90 J = 1, M - 1 
            DO 90 I = 1, N - 1 
   90    B (K, 3 * J, 3 * I) = (D (K, J - 1, I - 1) + 4.0D0 * D (K, J,  &
         I - 1) + D (K, J + 1, I - 1) + 4.0D0 * D (K, J - 1, I) +       &
         16.0D0 * D (K, J, I) + 4.0D0 * D (K, J + 1, I) + D (K, J - 1,  &
         I + 1) + 4.0D0 * D (K, J, I + 1) + D (K, J + 1, I + 1) )       &
         / 36.0D0                                                       
  100 END DO 
      RETURN 
      END SUBROUTINE BEZPNT                         
!                                                                       
!                                                                       
      SUBROUTINE BEZBRD (B, D, M, N) 
!*****************************************************************      
!                                                                *      
!     Auxiliary program for SUBROUTINE BEZIER.                   *      
!     The boundary BEZIER points are determined from the prede-  *      
!     termined boundary weight points using 1-dimensional BEZIER *      
!     splines.                                                   *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     M             INTEGER   number of patches in 1st direction *      
!     N             INTEGER   number ov patches in 2nd direction *      
!                                                                *      
!     D(3,0:M,0:N)  Double Precision coordinates of the          *      
!                   interpolation points                         *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B(3,0:3M,0:3N) Double Precision coordinates of the         *      
!                    boundary BEZIER points                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hartmut Turowski                                   *      
!  date     : 05.18.1990                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION B (3, 0:3 * M, 0:3 * N), D (3, 0:M, 0:N) 
      DO 50 K = 1, 3 
         DO 20 I = 0, N, N 
            DO 10 J = 1, M - 1 
               B (K, 3 * J - 2, 3 * I) = (2.0D0 * D (K, J - 1, I)       &
               + D (K, J, I) ) / 3.0D0                                  
               B (K, 3 * J, 3 * I) = (D (K, J - 1, I) + 4.0D0 * D (K, J,&
               I) + D (K, J + 1, I) ) / 6.0D0                           
               B (K, 3 * J + 2, 3 * I) = (D (K, J, I) + 2.0D0 * D (K, J &
               + 1, I) ) / 3.0D0                                        
   10       END DO 
            B (K, 2, 3 * I) = (D (K, 0, I) + 2.0D0 * D (K, 1, I) )      &
            / 3.0D0                                                     
            B (K, 3 * M - 2, 3 * I) = (2.0D0 * D (K, M - 1, I) + D (K,  &
            M, I) ) / 3.0D0                                             
   20    END DO 
         DO 40 J = 0, M, M 
            DO 30 I = 1, N - 1 
               B (K, 3 * J, 3 * I - 2) = (2.0D0 * D (K, J, I - 1)       &
               + D (K, J, I) ) / 3.0D0                                  
               B (K, 3 * J, 3 * I) = (D (K, J, I - 1) + 4.0D0 * D (K, J,&
               I) + D (K, J, I + 1) ) / 6.0D0                           
               B (K, 3 * J, 3 * I + 2) = (D (K, J, I) + 2.0D0 * D (K, J,&
               I + 1) ) / 3.0D0                                         
   30       END DO 
            B (K, 3 * J, 2) = (D (K, J, 0) + 2.0D0 * D (K, J, 1) )      &
            / 3.0D0                                                     
            B (K, 3 * J, 3 * N - 2) = (2.0D0 * D (K, J, N - 1) + D (K,  &
            J, N) ) / 3.0D0                                             
   40    END DO 
         B (K, 0, 0) = D (K, 0, 0) 
         B (K, 3 * M, 0) = D (K, M, 0) 
         B (K, 0, 3 * N) = D (K, 0, N) 
         B (K, 3 * M, 3 * N) = D (K, M, N) 
   50 END DO 
      RETURN 
      END SUBROUTINE BEZBRD                         

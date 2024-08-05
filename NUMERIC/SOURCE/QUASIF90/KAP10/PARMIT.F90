![  {Parametric Hermite Splines}                                        
![  {Parametric Hermite Splines}*)                                      
      SUBROUTINE PARMIT (N, MARG, X, Y, XDIREC, YDIREC, IDIREC, BCND1,  &
      BCNDN, AX, BX, CX, DX, EX, FX, AY, BY, CY, DY, EY, FY, T, IERR,   &
      XT, YT, SUP, DXT, AINF, PRC, AR1, AR2, AR3, H)                    
!                                                                       
!*****************************************************************      
!                                                                *      
!     PARMIT computes the coefficients of a parametric hermite   *      
!     spline. The end point conditions can be specified via MARG.*      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N       number of nodes  (X(I),Y(I))                       *      
!     MARG    index for the endpoint condition:                  *      
!               MARG = 1 : Periodic spline                       *      
!               MARG = 2 : Natural spline                        *      
!               MARG = 3 : User specified second derivatives at  *      
!                          the end points.                       *      
!                          In this case the program expects that *      
!                           2    2                               *      
!                          D Y/DX ( X(1) ) is in BCND1(1)  and   *      
!                           2    2                               *      
!                          D Y/DX ( X(N) )  in BCNDN(1).         *      
!               MARG = 4 : The user specifies the second         *      
!                          derivative                            *      
!                          ..     ..                             *      
!                          SX(T), SY(T) of the component splines *      
!                          at the end points.                    *      
!                          In this case the program expects to   *      
!                          find                                  *      
!                          ..                                    *      
!                          SX( T(1) )  in BCND1(1)               *      
!                          ..                                    *      
!                          SY( T(1) )  in BCND1(2)               *      
!                          ..                                    *      
!                          SX( T(N) )  in BCNDN(1), and          *      
!                          ..                                    *      
!                          SY( T(N) )  in BCNDN(2).              *      
!               MARG = 5 : The user specifies the curvature radii*      
!                          R1, RN at the end points.             *      
!                          In this case the program expects      *      
!                          R1  in BCND1(1), and                  *      
!                          RN  in BCNDN(1).                      *      
!                          Re. concavity : If the radius is      *      
!                          positive, the curvature circle lies   *      
!                          to the left of the spline viewed in   *      
!                          direction of increasing parameter     *      
!                          values (concave to the left); if the  *      
!                          radius is negative, the spline is     *      
!                          concave to the right.                 *      
!               MARG = 6 : The user specifies the third          *      
!                          derivative                            *      
!                          ...    ...                            *      
!                          SX(T), SY(T) of the component splines *      
!                          at the end points.                    *      
!                          In this case the program expects to   *      
!                          find                                  *      
!                          ...                                   *      
!                          SX( T(1) )  in BCND1(1)               *      
!                          ...                                   *      
!                          SY( T(1) )  in BCND1(2)               *      
!                          ...                                   *      
!                          SX( T(N) )  in BCNDN(1), and          *      
!                          ...                                   *      
!                          SY( T(N) )  in BCNDN(2).              *      
!     X       ) vectors ..(1:N); the given nodes (X(I),Y(I)) for *      
!     Y       ) I = 1, ... , N.                                  *      
!     XDIREC  vector XDIR(1:N);                                  *      
!             X components of the tangent or normal vector at    *      
!             ( X(I),Y(I) ), I=1,...,N, if the user prescribes   *      
!             such.                                              *      
!             Otherwise: not used                                *      
!     YDIREC  vector YDIR(1:N);                                  *      
!             Y components of the tangent or normal vector at    *      
!             ( X(I),Y(I) ), I=1,...,N, the value of DY/DX(X(I)) *      
!             for I =1,...,N; see parameter IDIREC. If one of the*      
!             values exceeds 1.0D+38, the program assumes a      *      
!             vertical tangent there.                            *      
!     IDIREC  index for specifying tangent/normal vectors:       *      
!               IDIREC = 1 : Tangent vectors given               *      
!               IDIREC = 2 : Normal vectors given                *      
!               IDIREC = 3 : derivatives DY/DX(X(I)) given       *      
!             Refer to XDIREC and YDIREC.                        *      
!             NOTE: When prescribing tangent or normal vectors   *      
!                   please note:                                 *      
!                   a) the length of such a vector does not in-  *      
!                      fluence the spline, since these vectors   *      
!                      are normalized internally.                *      
!                   b) the normed tangent vector will be stored  *      
!                      for all values of IDIREC in XT (X-        *      
!                      components) and YT (Y-components).        *      
!     BCND1   )  vectors ..(1:2);                                *      
!     BCNDN   )  end point conditions as given by the user via   *      
!             )  MARG.                                           *      
!             )  not used if not prescribed by the user.         *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     AX      )  vectors ..(1:N);                                *      
!     FX      )  the coefficients of the spline component SX(T)  *      
!             )                  .                               *      
!             )  for ( T(I),X(I),X(I) )                          *      
!                                                                *      
!     AY      )  vectors .. (1:N);                               *      
!     FY      )  the coefficients of the spline component SY(T)  *      
!             )                  .                               *      
!             )  for ( T(I),Y(I),Y(I) )                          *      
!                                                                *      
!     T       vector T(1:N); the parameter values T(I), I=1,...,N*      
!     IERR    error code.                                        *      
!               IERR = 0 : no error                              *      
!               IERR > 0 : error in input data, no output.       *      
!               Specifically:                                    *      
!               IERR = 1 : MARG < 1  or  MARG > 6                *      
!               IERR = 2 : IDIREC < 1  or  IDIREC > 3            *      
!               IERR = 3 : N < 3                                 *      
!               IERR = 4 : If IDIREC=3: For one point we cannot  *      
!                          compute a meaningful tangent due to   *      
!                          contradicting input. This will occur  *      
!                          e.g. if three adjacent points are col-*      
!                          linear with the central tangent speci-*      
!                          fied perpendicular to this line.      *      
!               IERR = 5 : Two cosecutive points coincide:       *      
!                          X(I)=X(I+1) and Y(I)=Y(I+1) for one   *      
!                          index I. (Non consecutive identical   *      
!                          points such as self crossing curves   *      
!                          are, however allowed.)                *      
!               IERR = 6 : a tangent or normal vector is equal to*      
!                          the zero vector.                      *      
!                           2     2                              *      
!               IERR = 7 : D  Y/DX   was specified at the end-   *      
!                          points, but the spline has a vertical *      
!                          tangent at one end point, leading to  *      
!                          mathematical inconsistencies.         *      
!               IERR = 8 : One of the curvature radii is zero.   *      
!               IERR = 9 : A periodic spline was stipulated, but *      
!                          the data satisfies  X(1).NE.X(N).     *      
!               IERR =10 : A periodic spline was stipulated, but *      
!                          the data satisfies  Y(1).NE.Y(N).     *      
!                                                                *      
!                                                                *      
!     AUXILIARY PARAMETERS:                                      *      
!     =====================                                      *      
!     XT      ) auxiliary vectors ..(1:N);                       *      
!     YT      ) the vectors XT and YT will contain the normalized*      
!     SUP     ) tangent vectors, see IDIREC.                     *      
!     DXT     )                                                  *      
!     AINF    )                                                  *      
!     PRC     )                                                  *      
!     AR1     )                                                  *      
!     AR2     )                                                  *      
!     AR3     )                                                  *      
!     H       )                                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: HERMIT                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Elmar Pohl                                         *      
!  Date     : 09.28.1985                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION X (N), Y (N), XDIREC (N), YDIREC (N), T (N),      &
      XT (N), YT (N)                                                    
      DOUBLEPRECISION BCND1 (2), BCNDN (2), SUP (N), DXT (N), AR1 (N),  &
      AR2 (N), AR3 (N)                                                  
      DOUBLEPRECISION AX (N), BX (N), CX (N), DX (N), EX (N), FX (N),   &
      H (N), PRC (N)                                                    
      DOUBLEPRECISION AY (N), BY (N), CY (N), DY (N), EY (N), FY (N),   &
      AINF (N)                                                          
      IERR = 0 
      IF (MARG.GE.1.AND.MARG.LE.6) GOTO 1 
      IERR = 1 
      RETURN 
    1 IF (IDIREC.GE.1.AND.IDIREC.LE.3) GOTO 2 
      IERR = 2 
      RETURN 
    2 IF (N.GE.3) GOTO 4 
      IERR = 3 
      RETURN 
    4 T (1) = 0.0D0 
      DO 6 I = 2, N 
         DELTX = X (I) - X (I - 1) 
         DELTY = Y (I) - Y (I - 1) 
         DELT = DELTX * DELTX + DELTY * DELTY 
         IF (DELT.GT.0.0D0) GOTO 5 
         IERR = 5 
         RETURN 
    5    T (I) = T (I - 1) + DSQRT (DELT) 
    6 END DO 
      IF (IDIREC - 2) 101, 102, 103 
  101 DO 19 I = 1, N 
         XT (I) = XDIREC (I) 
         YT (I) = YDIREC (I) 
   19 END DO 
      GOTO 104 
  102 DO 7 I = 1, N 
         XT (I) = YDIREC (I) 
         YT (I) = - XDIREC (I) 
    7 END DO 
      GOTO 104 
  103 DO 9 I = 1, N 
         IF (DABS (YDIREC (I) ) .GE.1.0D38) GOTO 8 
         XT (I) = 1.0D0 
         YT (I) = YDIREC (I) 
         GOTO 9 
    8    XT (I) = 0.0D0 
         YT (I) = 1.0D0 
    9 END DO 
  104 CONTINUE 
!                                                                       
      DO 11 I = 1, N 
         VLONG = DSQRT (XT (I) * XT (I) + YT (I) * YT (I) ) 
         IF (VLONG.GT.0.0D0) GOTO 10 
         IERR = 6 
         RETURN 
   10    XT (I) = XT (I) / VLONG 
         YT (I) = YT (I) / VLONG 
   11 END DO 
      NM1 = N - 1 
      IF ( (X (2) - X (1) ) * XT (1) + (Y (2) - Y (1) ) * YT (1) ) 25,  &
      27, 26                                                            
   27 IF (IDIREC - 3) 26, 28, 26 
   25 XT (1) = - XT (1) 
      YT (1) = - YT (1) 
   26 DO 12 I = 1, NM1 
         IF ( (X (I + 1) - X (I) ) * XT (I) + (Y (I + 1) - Y (I) )      &
         * YT (I) ) 20, 21, 12                                          
   21    IF ( (X (I) - X (I - 1) ) * XT (I) + (Y (I) - Y (I - 1) )      &
         * YT (I) ) 20, 29, 12                                          
   29    IF (IDIREC - 3) 12, 28, 12 
   20    XT (I) = - XT (I) 
         YT (I) = - YT (I) 
   12 END DO 
      IF ( (X (N) - X (NM1) ) * XT (N) + (Y (N) - Y (NM1) ) * YT (N) )  &
      22, 23, 24                                                        
   23 IF (IDIREC - 3) 24, 28, 24 
   28 IERR = 4 
      RETURN 
   22 XT (N) = - XT (N) 
      YT (N) = - YT (N) 
   24 GOTO (201, 202, 203, 204, 205, 206), MARG 
  201 MARGH = 1 
      GOTO 16 
  202 MARGH = 2 
      GOTO 16 
  203 IF (XT (1) .NE.0.0D0.AND.XT (N) .NE.0.0D0) GOTO 13 
      IERR = 7 
      RETURN 
   13 RB1 = 1.0D0 
      RBN = 1.0D0 
      GOTO 15 
  204 RB1 = BCND1 (1) 
      RBN = BCNDN (1) 
      GOTO 15 
  205 IF (BCND1 (1) .NE.0.0D0.AND.BCNDN (1) .NE.0.0D0) GOTO 14 
      IERR = 8 
      RETURN 
   14 RB1 = 1.0D0 
      IF (XT (1) .EQ.0.0D0) RB1 = - 1.0D0 / BCND1 (1) / YT (1) 
      RBN = 1.0D0 
      IF (XT (N) .EQ.0.0D0) RBN = - 1.0D0 / BCNDN (1) / YT (N) 
   15 MARGH = 3 
      GOTO 16 
  206 MARGH = 5 
      RB1 = BCND1 (1) 
      RBN = BCNDN (1) 
   16 IREP = 0 
      CALL HERMIT (N, MARGH, T, X, XT, RB1, RBN, IREP, AX, BX, CX, DX,  &
      EX, FX, MORSH, H, SUP, AINF, PRC, DXT, AR1, AR2, AR3)             
      IF (MARG.NE.1.OR.MORSH.NE.4) GOTO 17 
      IERR = 9 
      RETURN 
   17 GOTO (18, 18, 303, 304, 305, 306), MARG 
  303 RB1 = (XT (1) **3 * BCND1 (1) + YT (1) ) / XT (1) 
      RBN = (XT (N) **3 * BCNDN (1) + YT (N) ) / XT (N) 
      GOTO 18 
  304 RB1 = BCND1 (2) 
      RBN = BCNDN (2) 
      GOTO 18 
  305 RB1 = 1.0D0 
      IF (XT (1) .NE.0.0D0) RB1 = (1.0D0 / BCND1 (1) + YT (1) ) / XT (1) 
      RBN = 1.0D0 
      IF (XT (N) .NE.0.0D0) RBN = (1.0D0 / BCNDN (1) + YT (N) ) / XT (N) 
      GOTO 18 
  306 RB1 = BCND1 (2) 
      RBN = BCNDN (2) 
   18 IREP = 1 
      CALL HERMIT (N, MARGH, T, Y, YT, RB1, RBN, IREP, AY, BY, CY, DY,  &
      EY, FY, MORSH, H, SUP, AINF, PRC, DXT, AR1, AR2, AR3)             
      IF (MARG.EQ.1.AND.MORSH.EQ.4) IERR = 10 
      RETURN 
      END SUBROUTINE PARMIT                         
!                                                                       
!                                                                       
      SUBROUTINE PMTVAL (N, T0, T, AX, BX, CX, DX, EX, FX, AY, BY, CY,  &
      DY, EY, FY, SX, SY, OUTP)                                         
!                                                                       
!*******************************************************************    
!                                                                  *    
!     PMTVAL computes functional values of a parametric hermite    *    
!     polynomial spline of fifth degree ( SX(T),SY(T) ) and its    *    
!     derivatives at T=T0.                                         *    
!     While this program could be used to obtain an equidistant    *    
!     table of values for the spline for graphing e.g., this is not*    
!     recommended, since PMTVAL performs an expensive interval     *    
!     search for each input T0. Moreover not all derivatives will  *    
!     generally be required.                                       *    
!     To make a table of values we recommend a program like PMTAB. *    
!                                                                  *    
!                                                                  *    
!     INPUT PARAMETERS:                                            *    
!     =================                                            *    
!     T0      value where we want to evaluate the spline           *    
!     N       number of parameter values T(I)                      *    
!     T       vector T(1:N); parameter values of the splines,      *    
!             as computed by the SUBROUTINE PARMIT for example.    *    
!     AX      )                                                    *    
!     BX      )  vectors ..(1:N);                                  *    
!     CX      )  the coefficients of the spline component SX       *    
!     DX      )                                                    *    
!     EX      )                                                    *    
!     FX      )                                                    *    
!                                                                  *    
!     AY      )                                                    *    
!     BY      )  vectors .. (1:N);                                 *    
!     CY      )  the coefficients of the spline component SX       *    
!     DY      )                                                    *    
!     EY      )                                                    *    
!     FY      )                                                    *    
!                                                                  *    
!                                                                  *    
!     OUTPUT PARAMETERS:                                           *    
!     ==================                                           *    
!     SX      SX = SX(T0)                                          *    
!     SY      SY = SY(T0)                                          *    
!     OUTP    2-dimensional array OUTP(1:5,1:2) for the derivatives*    
!                           (K)                                    *    
!             OUTP(K,1) = SX    (T0)                               *    
!                           (K)                                    *    
!             OUTP(K,2) = SY    (T0),  K=1(1)5                     *    
!                                                                  *    
!------------------------------------------------------------------*    
!                                                                  *    
!  Subroutines required: HMTVAL                                    *    
!                                                                  *    
!*******************************************************************    
!                                                                  *    
!  Author   : Elmar Pohl                                           *    
!  Date     : 28.09.1985                                           *    
!  Source   : FORTRAN 77                                           *    
!                                                                  *    
!*******************************************************************    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION T (N), AX (N), BX (N), CX (N), DX (N), EX (N),    &
      FX (N), AY (N)                                                    
      DOUBLEPRECISION BY (N), CY (N), DY (N), EY (N), FY (N), AUSG (5), &
      OUTP (5, 2)                                                       
      DOUBLEPRECISION HMTVAL 
      SX = HMTVAL (N, T0, AX, BX, CX, DX, EX, FX, T, AUSG) 
      DO 10 I = 1, 5 
         OUTP (I, 1) = AUSG (I) 
   10 END DO 
      SY = HMTVAL (N, T0, AY, BY, CY, DY, EY, FY, T, AUSG) 
      DO 20 I = 1, 5 
         OUTP (I, 2) = AUSG (I) 
   20 END DO 
      RETURN 
      END SUBROUTINE PMTVAL                         

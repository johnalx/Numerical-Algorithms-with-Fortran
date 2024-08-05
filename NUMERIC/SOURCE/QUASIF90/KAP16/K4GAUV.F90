      SUBROUTINE K4GAUV (USERF, X, NX, Y, NY, METHOD, MOLD, CREC,       &
      ESTDIV, DIVIAT, WORK, IERR, IUFCLL)                               
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over rectangles using NEWTON-COTES formulas for       *      
! GAUSSIAN nodes:                                                *      
!                                                                *      
! The FUNCTION USERF(X,Y) shall be integrated using summed       *      
! GAUSSIAN formulas for the rectangle [A,B] x [C,D].             *      
! The sub-rectangles to be used are given via the vectors        *      
! X and Y.                                                       *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! USERF   : user defined FUNCTION USERF(X,Y), whose integral is  *      
!           to be computed.                                      *      
!           The FUNCTION USERF must be declared as EXTERNAL in   *      
!           the calling program.                                 *      
!           The FUNCTION should have the following form:         *      
!                  DOUBLE PRECISION FUNCTION USERF(X,Y)          *      
!                  DOUBLE PRECISION X,Y                          *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  USERF=F(X,Y)                                  *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  RETURN                                        *      
!                  STOP                                          *      
! X       : DOUBLE PRECISION vector X(0:NX) containing the X     *      
!           partition of [A,B].                                  *      
!           A = X(0) < X(1) < ... < X(NX) = B                    *      
! NX      : number of intervals in X-direction, NX > 0           *      
! Y       : DOUBLE PRECISION vector Y(0:NY) containing the Y     *      
!           partition of [C,D].                                  *      
!           C = Y(0) < Y(1) < ... < Y(NY) = D                    *      
! NY      : number of intervals in Y-direction, NY > 0           *      
! METHOD  : INTEGER indicating method used, 0 <= METHOD <= 7     *      
! MOLD    : INTEGER, the number in METHOD at the previous call.  *      
!           Upon first call we must have: MOLD different from    *      
!           METHOD                                               *      
!           If K4NECN is called repeatedly with METHOD=MOLD the  *      
!           internal initializing of parameters is skipped.      *      
! ESTDIV  : LOGICAL variable, indicates whether error estimate   *      
!           is to be computed (ESTDIV=TRUE) or not (ESTDIV=FALSE)*      
! WORK    : 2-dimensional DOUBLE PRECISION array                 *      
!           WORK(3,0:METHOD-1). If METHOD=MOLD this array must   *      
!           contain the initializing parameters for the method.  *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! MOLD    : INTEGER indicating method used                       *      
! CREC    : DOUBLE PRECISION value for the integral              *      
! WORK    : 2-dimensional DOUBLE PRECISION array WORK(2,0:METHOD)*      
!           containing the constants for the method              *      
! DIVIAT  : DOUBLE PRECISION error estimate                      *      
!           If ESTDIV=TRUE the error is estimated by one extra   *      
!           cubature for the halved step size.                   *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 X-interval of length zero    *      
!                            IERR=2 Y-interval of length zero    *      
!                            IERR=3 Number of method erroneous   *      
!                            IERR=4 NX < 1 or NY < 1             *      
! IUFCLL  : INTEGER, the number of functional evaluations        *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! =================                                              *      
! I, J, K : INTEGER, loop counters                               *      
! II, JJ  : INTEGER, loop counters                               *      
! I1, J1  : INTEGER, loop counters                               *      
! KMAX    : INTEGER number of cubature passes                    *      
! DBLEX   : DOUBLE PRECISION value for K                         *      
! HX      : DOUBLE PRECISION step size in X-direction            *      
! HY      : DOUBLE PRECISION step size in Y-direction            *      
! HXM     : DOUBLE PRECISION mid-point of X-interval             *      
! HYM     : DOUBLE PRECISION mid-point of Y-interval             *      
! CRECH   : DOUBLE PRECISION variable used for error estimation  *      
! FAC     : DOUBLE PRECISION variable used for CREC              *      
! HELPF   : DOUBLE PRECISION variable used for CREC              *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: K4GINI                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 06.12.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (2, 0:METHOD), X (0:NX), Y (0:NY), CREC,     &
      CRECH, DIVIAT, HX, HY, HXM, HYM, DBLEK, FAC, HELPF, USERF         
!                                                                       
! LOGICAL variable ESTDIV                                               
!                                                                       
      LOGICAL ESTDIV 
!                                                                       
! Initialize IUFCLL                                                     
!                                                                       
      IUFCLL = 0 
!                                                                       
! Validate input data                                                   
!                                                                       
!   Length of X-intervals                                               
!                                                                       
      DO 10 I = 1, NX 
         IF (X (I) .LE.X (I - 1) ) THEN 
            IERR = 1 
            RETURN 
         ENDIF 
   10 END DO 
!                                                                       
!   Length of Y-intervals                                               
!                                                                       
      DO 20 I = 1, NY 
         IF (Y (I) .LE.Y (I - 1) ) THEN 
            IERR = 2 
            RETURN 
         ENDIF 
   20 END DO 
!                                                                       
!   Check number of method                                              
!                                                                       
      IF (METHOD.LT.0.OR.METHOD.GT.7) THEN 
         IERR = 3 
         RETURN 
!                                                                       
!   Check number of sub-intervals                                       
!                                                                       
      ELSEIF (NX.LE.0.OR.NY.LE.0) THEN 
         IERR = 4 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
! If necessary, check initial values                                    
!                                                                       
      IF (METHOD.NE.MOLD) THEN 
         CALL K4GINI (METHOD, WORK) 
         MOLD = METHOD 
      ENDIF 
!                                                                       
!                                                                       
      IF (ESTDIV) THEN 
         KMAX = 2 
      ELSE 
         KMAX = 1 
      ENDIF 
!                                                                       
! Loop over necessary cubature runs                                     
!                                                                       
      DO 30 K = 1, KMAX 
!                                                                       
! Change K                                                              
!                                                                       
         DBLEK = DBLE (K) 
!                                                                       
! Initialize CREC                                                       
!                                                                       
         CREC = 0.0D0 
!                                                                       
! Find approximation for the integral                                   
!                                                                       
         DO 40 I = 0, NX - 1 
!                                                                       
! Find step size in X-direction                                         
!                                                                       
            HX = (X (I + 1) - X (I) ) / (2.0D0 * DBLEK) 
            DO 50 I1 = 1, 2 * K - 1, 2 
               HXM = X (I) + DBLE (I1) * HX 
               DO 60 J = 0, NY - 1 
!                                                                       
! Find step size in Y-direction                                         
!                                                                       
                  HY = (Y (I + 1) - Y (I) ) / (2.0D0 * DBLEK) 
                  DO 70 J1 = 1, 2 * K - 1, 2 
                     HYM = Y (J) + DBLE (J1) * HY 
                     DO 80 II = 0, METHOD 
                        DO 90 JJ = 0, METHOD 
                           FAC = HX * HY * WORK (2, II) * WORK (2, JJ) 
                           HELPF = USERF (HXM + HX * WORK (1, II),      &
                           HYM + HY * WORK (1, JJ) )                    
                           IUFCLL = IUFCLL + 1 
                           CREC = CREC + FAC * HELPF 
   90                   END DO 
   80                END DO 
   70             END DO 
   60          END DO 
   50       END DO 
   40    END DO 
!                                                                       
! When estimating the error, store first integral value                 
!                                                                       
         IF (ESTDIV.AND.K.EQ.1) CRECH = CREC 
   30 END DO 
!                                                                       
! Error estimation                                                      
!                                                                       
      IF (ESTDIV) DIVIAT = (CREC - CRECH) / 3.0D0 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4GAUV                         

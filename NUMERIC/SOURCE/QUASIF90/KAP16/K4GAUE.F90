![          {Gau"s Cubature Formulas for Rectangles}*)                  
      SUBROUTINE K4GAUE (USERF, A, B, IP, C, D, IQ, METHOD, MOLD, CREC, &
      ESTDIV, DIVIAT, WORK, IERR, IUFCLL)                               
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over rectangles via NEWTON-COTES formulas for GAUSSIAN*      
! nodes.                                                         *      
!                                                                *      
! The FUNCTION USERF(X,Y) shall be integrated using summed       *      
! GAUSSIAN formulas for the rectangle [A,B] x [C,D].             *      
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
!                                                                *      
! A       : DOUBLE PRECISION left hand endpoint in X-direction   *      
! B       : DOUBLE PRECISION right hand endpoint in X-direction  *      
! IP      : INTEGER, the number of sub-intervals in X-direction  *      
! C       : DOUBLE PRECISION lower endpoint in Y-direction       *      
! D       : DOUBLE PRECISION upper endpoint in Y-direction       *      
! IQ      : INTEGER, the number of sub-intervals in Y-direction  *      
! METHOD  : INTEGER, designating the method: 0 <= METHOD <= 7 :  *      
!           METHOD = 0 : trapezoidal rule                        *      
!                  = 1 : summed trapezoidal rule                 *      
!                  = 2 : Simpson's rule                          *      
!                  = 3 : 3/8 rule                                *      
!                  = 4 : 4/90 rule                               *      
!                  = 5 : 5/288 rule                              *      
!                  = 6 : 6/840 rule                              *      
!                  = 7 : 7/17280 rule.                           *      
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
! DIVIAT  : DOUBLE PRECISION error estimate                      *      
!           If ESTDIV=TRUE the error is estimated by one extra   *      
!           cubature for the halved step size.                   *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 number of intervals in       *      
!                                   X-direction erroneous        *      
!                            IERR=2 number of intervals in       *      
!                                   Y-direction erroneous        *      
!                            IERR=3 Number of method erroneous   *      
!                            IERR=4 interval of integration has  *      
!                                   length zero                  *      
! IUFCLL  : INTEGER, the number of functional evaluations        *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! =================                                              *      
! I, J, K : INTEGER, loop counters                               *      
! II, JJ  : INTEGER, loop counters                               *      
! KMAX    : INTEGER number of cubature passes                    *      
! IPX     : INTEGER number of intervals in X-direction           *      
! IPY     : INTEGER number of intervals in Y-direction           *      
! DI      : DOUBLE PRECISION value for 2*I+1                     *      
! DJ      : DOUBLE PRECISION value for 2*J+1                     *      
! HX      : DOUBLE PRECISION step size in X-direction            *      
! HY      : DOUBLE PRECISION step size in Y-direction            *      
! CRECH   : DOUBLE PRECISION variable used for error estimation  *      
! FAC     : DOUBLE PRECISION variable used for CREC              *      
! HELPC   : DOUBLE PRECISION variable used for CREC              *      
! HELPX   : DOUBLE PRECISION variable used for CREC              *      
! HELPY   : DOUBLE PRECISION variable used for CREC              *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: K4GINI                                  *      
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
      DOUBLEPRECISION WORK (2, 0:METHOD), A, B, C, D, CREC, CRECH,      &
      DIVIAT, HX, HY, DI, DJ, FAC, HELPX, HELPY, HELPC, USERF           
!                                                                       
! Initialize the LOGICAL variable ESTDIV                                
!                                                                       
      LOGICAL ESTDIV 
!                                                                       
! Initialize IUFCLL                                                     
!                                                                       
      IUFCLL = 0 
!                                                                       
! Check input data for validity                                         
!                                                                       
      IF (IP.LT.1) THEN 
         IERR = 1 
         RETURN 
      ELSEIF (IQ.LT.1) THEN 
         IERR = 2 
         RETURN 
      ELSEIF (METHOD.LT.0.OR.METHOD.GT.7) THEN 
         IERR = 3 
         RETURN 
      ELSEIF (A.EQ.B.OR.C.EQ.D) THEN 
         IERR = 4 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
! Initialize as necessary                                               
!                                                                       
      IF (METHOD.NE.MOLD) THEN 
         CALL K4GINI (METHOD, WORK) 
         MOLD = METHOD 
      ENDIF 
!                                                                       
! Determine number of needed cubature passes                            
!                                                                       
      IF (ESTDIV) THEN 
         KMAX = 2 
      ELSE 
         KMAX = 1 
      ENDIF 
!                                                                       
! Determine actual number of sub-intervals                              
! in X- and Y-directions                                                
!                                                                       
      DO 10 K = 1, KMAX 
         IPX = K * IP 
         IQY = K * IQ 
!                                                                       
! Initialize CREC                                                       
!                                                                       
         CREC = 0.0D0 
!                                                                       
! Determine step sizes in X- and Y-directions                           
!                                                                       
         HX = 0.5D0 * (B - A) / DBLE (IPX) 
         HY = 0.5D0 * (D-C) / DBLE (IQY) 
!                                                                       
! Find approximation for the integral                                   
!                                                                       
         DO 20 I = 0, IPX - 1 
            DI = 2.0D0 * DBLE (I) + 1.0D0 
            DO 30 J = 0, IQY - 1 
               DJ = 2.0D0 * DBLE (J) + 1.0D0 
               DO 40 II = 0, METHOD 
                  DO 50 JJ = 0, METHOD 
                     FAC = HX * HY * WORK (2, II) * WORK (2, JJ) 
                     HELPX = A + HX * (WORK (1, II) + DI) 
                     HELPY = C + HY * (WORK (1, JJ) + DJ) 
                     HELPC = USERF (HELPX, HELPY) 
                     IUFCLL = IUFCLL + 1 
                     CREC = CREC + FAC * HELPC 
   50             END DO 
   40          END DO 
   30       END DO 
   20    END DO 
!                                                                       
! If estimating the error, store the first integral value               
!                                                                       
         IF (ESTDIV.AND.K.EQ.1) CRECH = CREC 
   10 END DO 
!                                                                       
! Error estimation                                                      
!                                                                       
      IF (ESTDIV) DIVIAT = (CREC - CRECH) / 3.0D0 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4GAUE                         
!                                                                       
!                                                                       
      SUBROUTINE K4GINI (METHOD, WORK) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Subroutine that initializes the constants in K4GAUE depending  *      
! on method chosen.                                              *      
!                                                                *      
!                                                                *      
! INPUT PARAMETER:                                               *      
! ================                                               *      
! METHOD  : INTEGER designating the method, 0 <= METHOD <= 7     *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! WORK    : 2-dimensional DOUBLE PRECISION array                 *      
!           WORK(2,0:METHOD) containing the constants for the    *      
!           method                                               *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLE:                                                *      
! ===============                                                *      
! I       : INTEGER, loop parameter                              *      
! J       : INTEGER, auxiliary variable                          *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 12.06.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (2, 0:METHOD) 
!                                                                       
! Set up upper half of WORK                                             
!                                                                       
      IF (METHOD.EQ.0) THEN 
         WORK (1, 0) = 0.0D0 
         WORK (2, 0) = 2.0D0 
         J = 0 
      ELSEIF (METHOD.EQ.1) THEN 
         WORK (1, 0) = - 0.577350269189626D0 
         WORK (2, 0) = 1.0D0 
         J = 0 
      ELSEIF (METHOD.EQ.2) THEN 
         WORK (1, 0) = - 0.774596669241483D0 
         WORK (2, 0) = 0.5555555555555556D0 
         WORK (1, 1) = 0.0D0 
         WORK (2, 1) = 0.8888888888888888D0 
         J = 1 
      ELSEIF (METHOD.EQ.3) THEN 
         WORK (1, 0) = - 0.861136311594053D0 
         WORK (2, 0) = 0.347854845137454D0 
         WORK (1, 1) = - 0.339981043584856D0 
         WORK (2, 1) = 0.652145154862546D0 
         J = 1 
      ELSEIF (METHOD.EQ.4) THEN 
         WORK (1, 0) = - 0.906179845938664D0 
         WORK (2, 0) = 0.236926885056189D0 
         WORK (1, 1) = - 0.538469310105683D0 
         WORK (2, 1) = 0.478628670499366D0 
         WORK (1, 2) = 0.0D0 
         WORK (2, 2) = 0.5688888888888889D0 
         J = 2 
      ELSEIF (METHOD.EQ.5) THEN 
         WORK (1, 0) = - 0.9324695142031521D0 
         WORK (2, 0) = 0.17132449237917D0 
         WORK (1, 1) = - 0.661209386466265D0 
         WORK (2, 1) = 0.360761573048139D0 
         WORK (1, 2) = - 0.238619186083197D0 
         WORK (2, 2) = 0.467913934572691D0 
         J = 2 
      ELSEIF (METHOD.EQ.6) THEN 
         WORK (1, 0) = - 0.949107912342759D0 
         WORK (2, 0) = 0.12948496616887D0 
         WORK (1, 1) = - 0.741531185599394D0 
         WORK (2, 1) = 0.279705391489277D0 
         WORK (1, 2) = - 0.405845151377397D0 
         WORK (2, 2) = 0.381830050505119D0 
         WORK (1, 3) = 0.0D0 
         WORK (2, 3) = 0.417959183673469D0 
         J = 3 
      ELSEIF (METHOD.EQ.7) THEN 
         WORK (1, 0) = - 0.960289856497536D0 
         WORK (2, 0) = 0.101228536290376D0 
         WORK (1, 1) = - 0.7966664774136269D0 
         WORK (2, 1) = 0.222381034453374D0 
         WORK (1, 2) = - 0.525532409916329D0 
         WORK (2, 2) = 0.313706645877887D0 
         WORK (1, 3) = - 0.18343464249565D0 
         WORK (2, 3) = 0.362683783378362D0 
         J = 3 
      ENDIF 
!                                                                       
! Set up lower half of WORK by symmetry                                 
!                                                                       
      DO 10 I = 0, J 
         WORK (1, METHOD-I) = - WORK (1, I) 
         WORK (2, METHOD-I) = WORK (2, I) 
   10 END DO 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4GINI                         

      SUBROUTINE HSTART (DES, N, X, BETA, Y, RELERR, ABSERR, QG, DSMALL,&
      DLARGE, H)                                                        
!                                                                       
!*****************************************************************      
!                                                                *      
! HSTART computes the initial step size for solving an initial   *      
! value problem numerically. The number of differential          *      
! equations in the system is limited to 12.                      *      
! To compute this step size we determine a LIPSCHITZ constant,   *      
! an upper bound for the first and second derivative of the      *      
! differential equation in a neighborhood of X=X0.               *      
! The algorithm used in HSTART is adapted from the software      *      
! package DEPAC (design of a user oriented package of ode        *      
! solvers).                                                      *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X       : DOUBLE PRECISION initial value for the integration:  *      
!           X=X0                                                 *      
! BETA    : DOUBLE PRECISION endpoint X=BETA at which we want to *      
!           find the solution                                    *      
! N       : number of differential equations in the system,      *      
!           or the size of Y:   0 < N < 13                       *      
! DES     : SUBROUTINE DES must be declared as EXTERNAL in the   *      
!           calling program. DES describes the system of         *      
!           differential equations and must have the following   *      
!           form:                                                *      
!                  SUBROUTINE DES(X,Y,N,YPUNKT)                  *      
!                  DOUBLE PRECISION Y(N),YPUNKT(N),X             *      
!                  YPUNKT(1)=....                                *      
!                  YPUNKT(2)=....                                *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                   YPUNKT(N)=....                               *      
!                   RETURN                                       *      
!                   END                                          *      
! Y       : DOUBLE PRECISION vector Y(1:N), the solution at X=X0 *      
! ABSERR  : DOUBLE PRECISION error bound for the absolute error  *      
!           (ABSERR >= 0). If ABSERR=0, then only the relative   *      
!           error is checked.                                    *      
! RELERR  : DOUBLE PRECISION error bound for the relative error  *      
!           (RELERR >= 0). If RELERR=0, then only the absolute   *      
!           error is checked.                                    *      
! QG      : global error order of the method                     *      
! DSMALL  : DOUBLE PRECISION machine constant                    *      
! DLARGE  : DOUBLE PRECISION largest representable number        *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! H       : DOUBLE PRECISION computed step size                  *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! WORK    : 2-dim. DOUBLE PRECISION array WORK(1:12,1:5)         *      
! DX      : DOUBLE PRECISION interval length                     *      
! ABSDX   : DOUBLE PRECISION absolute value of DX                *      
! RELPER  : DOUBLE PRECISION, RELPER=DSMALL**0.375D0           *        
! DA      : DOUBLE PRECISION variation in X                      *      
! DELF    : DOUBLE PRECISION auxiliary variable                  *      
! DFDXB   : DOUBLE PRECISION upper bound for second derivative,  *      
!           determined via differential quotient                 *      
! FBND    : DOUBLE PRECISION upper bound for first derivative    *      
! DELY    : DOUBLE PRECISION auxiliary variable                  *      
! DFDUB   : DOUBLE PRECISION LIPSCHITZ constant                  *      
! DY      : DOUBLE PRECISION auxiliary variable                  *      
! YDPB    : DOUBLE PRECISION upper bound of second derivative    *      
! TOLMIN  :-                                                     *      
! TOLSUM  :- DOUBLE PRECISION auxiliary variables for            *      
! TOL     :-                  computing TOLP                     *      
! TOLEXP  :-                                                     *      
! TOLP    : DOUBLE PRECISION tolerance value                     *      
! SRYDPB  : DOUBLE PRECISION, SRYDPB=SQRT(0.5*YDPB)              *      
! J, K    : loop variables                                       *      
! LK      : Number of iterations to compute the LIPSCHITZ        *      
!           constant                                             *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: VMNORM                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Volker KrÅger                                     *      
!  Date      : 07.08.1990                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION Y (N), WORK (12, 5) 
      DOUBLEPRECISION X, BETA, RELERR, ABSERR, QG, DSMALL, DLARGE, H,   &
      DX, ABSDX, RELPER, DA, DELF, VMNORM, DFDXB, FBND, DELY, DFDUB, DY,&
      YDPB, TOLMIN, TOLSUM, TOL, TOLEXP, TOLP, SRYDPB                   
!                                                                       
! Determine an upper bound for the second derivative (DFDXB) via        
! the differential quotient and an upper bound for the first            
! derivative (FBND).                                                    
!                                                                       
      DX = BETA - X 
      ABSDX = DABS (DX) 
      RELPER = DSMALL**0.375D0 
      DA = DSIGN (DMAX1 (DMIN1 (RELPER * DABS (X), ABSDX), 1.0D2 *      &
      DSMALL * DABS (X) ), DX)                                          
      IF (DABS (DA) .LT.DSMALL) DA = RELPER * DX 
      CALL DES (X + DA, Y, N, WORK (1, 1) ) 
      CALL DES (X, Y, N, WORK (1, 5) ) 
      DO 10 J = 1, N 
         WORK (J, 2) = WORK (J, 1) - WORK (J, 5) 
   10 END DO 
      DELF = VMNORM (WORK (1, 2), N) 
      DFDXB = DLARGE 
      IF (DELF.LT.DLARGE * DABS (DA) ) DFDXB = DELF / DABS (DA) 
      FBND = VMNORM (WORK (1, 1), N) 
!                                                                       
! Estimate the LIPSCHITZ constant (DFDUB) of the system of              
! differential equations and chose an upper bound (FBND)                
! for the first derivative.                                             
!                                                                       
      DELY = RELPER * VMNORM (Y, N) 
      IF (DELY.LT.DSMALL) DELY = RELPER 
      DELY = DSIGN (DELY, DX) 
      DELF = VMNORM (WORK (1, 5), N) 
      FBND = DMAX1 (FBND, DELF) 
      IF (DELF.LT.DSMALL) THEN 
         DO 40 J = 1, N 
            WORK (J, 3) = 0.0D0 
            WORK (J, 2) = 1.0D0 
   40    END DO 
         DELF = 1.0D0 
      ELSE 
         DO 20 J = 1, N 
            WORK (J, 3) = WORK (J, 5) 
            WORK (J, 2) = WORK (J, 5) 
   20    END DO 
      ENDIF 
      DFDUB = 0.0D0 
      LK = MIN (N + 1, 3) 
      DO 140 K = 1, LK 
         DO 60 J = 1, N 
            WORK (J, 4) = Y (J) + DELY * (WORK (J, 2) / DELF) 
   60    END DO 
         IF (K.EQ.2) THEN 
            CALL DES (X + DA, WORK (1, 4), N, WORK (1, 2) ) 
            DO 90 J = 1, N 
               WORK (J, 4) = WORK (J, 2) - WORK (J, 1) 
   90       END DO 
         ELSE 
            CALL DES (X, WORK (1, 4), N, WORK (1, 2) ) 
            DO 70 J = 1, N 
               WORK (J, 4) = WORK (J, 2) - WORK (J, 5) 
   70       END DO 
         ENDIF 
         FBND = DMAX1 (FBND, VMNORM (WORK (1, 2), N) ) 
         DELF = VMNORM (WORK (1, 4), N) 
         IF (DELF.GE.DLARGE * DABS (DELY) ) THEN 
            DFDUB = DLARGE 
            GOTO 150 
         ENDIF 
         DFDUB = DMAX1 (DFDUB, DELF / DABS (DELY) ) 
         IF (K.LT.LK) THEN 
            IF (DELF.LT.DSMALL) DELF = 1.0D0 
            DO 130 J = 1, N 
               IF (K.EQ.2) THEN 
                  DY = Y (J) 
                  IF (DABS (DY) .LT.DSMALL) DY = DELY / RELPER 
               ELSE 
                  DY = DABS (WORK (J, 4) ) 
                  IF (DY.LT.DSMALL) DY = DELF 
               ENDIF 
               IF (DABS (WORK (J, 3) ) .LT.DSMALL) WORK (J, 3) = WORK ( &
               J, 2)                                                    
               DY = DSIGN (DY, WORK (J, 3) ) 
               WORK (J, 2) = DY 
  130       END DO 
            DELF = VMNORM (WORK (1, 2), N) 
         ENDIF 
  140 END DO 
  150 YDPB = DFDXB + DFDUB * FBND 
!                                                                       
! Set tolerance value (TOLP) for computing the initial                  
! step size.                                                            
!                                                                       
      TOLMIN = DLARGE 
      TOLSUM = 0.0D0 
      DO 170 K = 1, N 
         TOL = RELERR * DABS (Y (K) ) + ABSERR 
         IF (TOL.LT.DSMALL) TOL = DABS (DELY) * RELERR 
         TOLEXP = DLOG10 (TOL) 
         TOLMIN = DMIN1 (TOLMIN, TOLEXP) 
         TOLSUM = TOLSUM + TOLEXP 
  170 END DO 
      TOLP = 1.0D1** (0.5D0 * (TOLSUM / DBLE (N) + TOLMIN) / (QG +      &
      1.0D0) )                                                          
!                                                                       
! Determine initial step size and direction of the integration          
!                                                                       
      H = ABSDX 
      IF (YDPB.GT.DSMALL.OR.FBND.GT.DSMALL) THEN 
         IF (YDPB.GT.DSMALL) THEN 
            SRYDPB = DSQRT (0.5D0 * YDPB) 
            IF (TOLP.LT.SRYDPB * ABSDX) H = TOLP / SRYDPB 
         ELSEIF (TOLP.LT.FBND * ABSDX) THEN 
            H = TOLP / FBND 
         ENDIF 
      ELSEIF (TOLP.LT.1.0D0) THEN 
         H = ABSDX * TOLP 
      ENDIF 
      IF (H * DFDUB.GT.1.0D0) H = 1.0D0 / DFDUB 
      H = DMAX1 (H, 1.0D2 * DSMALL * DABS (X) ) 
      IF (H.LT.DSMALL) H = DSMALL * DABS (BETA) 
      H = DSIGN (H, DX) 
      RETURN 
      END SUBROUTINE HSTART                         

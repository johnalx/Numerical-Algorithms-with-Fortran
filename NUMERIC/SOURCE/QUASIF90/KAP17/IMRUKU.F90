!                                                                       
!*****************************************************************      
!                                                                *      
! IMRUKU solves an initial value problem (IVP) for systems of N  *      
! given differential equations of 1st order.                     *      
! The solution of the IVP is obtained with an implicit RUNGE-    *      
! KUTTA method (IRKM), where we use step size control and control*      
! the order of the IRKM as well.                                 *      
! For this the subroutine IMRUKU needs IRKMs of orders 1 to MMAX.*      
! The coefficients that define these implicit RUNGE-KUTTA methods*      
! are assumed to be available in a file available at the input   *      
! LUNIN. This file can be created using the SUBROUTINE IRKCOE.   *      
!                                                                *      
! Before each RUNGE-KUTTA step we precheck the efficiency AW in  *      
! order to find the optimal order in relation to the number of   *      
! functional evaluations according to                            *      
!         AW(EPS,M) = (N+1+4*M**2)/H(EPS,M).                     *      
! Here H(EPS,M) is the step size for the order M and the accuracy*      
! bound EPS. For each step we chose an M for which the work AW   *      
! required for the IRKMs of order 1 to M+1 satisfy:              *      
!        AW(EPS,I) >  AW(EPS,I+1) for I=1, ..., M-1              *      
!  and   AW(EPS,M) <= AW(EPS,M+1);                               *      
!  or    M=MMAX-1, if there is no suitable M between 1 and MMAX-1*      
!                                                                *      
! For this M we chose a step size H which usually is taken equal *      
! to the theoretical step size H(EPS,M).                         *      
! The RUNGE-KUTTA step is performed for two IRKMs. One IRKM has  *      
! the order M with coefficient A, ALPHA and BETA. The other IRKM *      
! is of order M+1 for the coefficients AQ, ALPHAQ and BETAQ.     *      
! Determining the KI, called DBLEK and DBLEKQ in the program, is *      
! done iteratively.                                              *      
!                                                                *      
! If in an iteration step, the relative difference EREL of the   *      
! two approximations Y and YQ, that were obtained by two         *      
! different methods, exceeds EPS, then the step size H is reduced*      
! to H * SF * (0.5 * EPS * E)**(1/(2*M+1)). Here SF is a safety  *      
! factor and E is the absolute difference of Y and YQ. Then the  *      
! last step is repeated with the new step size H.                *      
! If in an iteration step, the difference between two successive *      
! approximations YQ(ALT) and YQ(NEU) exceeds the estimated       *      
! theoretical iteration error DK, no convergence can be expected.*      
! In this case the step size is decreased to 6/10 of the old     *      
! step size and the last integration step is repeated.           *      
! Iteration is continued until the relative difference between   *      
! two successive approximations becomes smaller than the required*      
! accuracy EPS.                                                  *      
! Based on theoretical considerations this should be the case    *      
! after at most 2*M+1 iterations.                                *      
! However, if this does not happen, the order is increased to    *      
! M + 1 and the last step is recomputed. If this would exceed    *      
! MMAX-1 we alter the step size to 0.8*H instead.                *      
! If, due to our step size control, H becomes less than ten      *      
! times the machine constant, we stop : no convergence can be    *      
! expected on the computer used.                                 *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! DES     : name of a subroutine, that represents the system of  *      
!           differential equations. This subroutine has to be    *      
!           declared as EXTERNAL in the calling program. It has  *      
!           the form:                                            *      
!                 SUBROUTINE DES (N,X,Y,F)                       *      
!                 INTEGER N                                      *      
!                 DOUBLE PRECISION X, Y(N), F(N)                 *      
!                 F(1) = F1(X,Y(1),Y(2),...,Y(N))                *      
!                 F(2) = F2(X,Y(1),Y(2),...,Y(N))                *      
!                 -------------------------------                *      
!                 F(N) = FN(X,Y(1),Y(2),...,Y(N))                *      
!                 RETURN                                         *      
!                 END                                            *      
! N       : The number of differential equations                 *      
! MMAX    : maximally allowed order, i.e. the highest order for  *      
!           which coefficients are available in an auxiliary     *      
!           file addressed at LUNIN. MMAX has to be >= 5.        *      
!           It is not advisable to chose the maximum order NMAX  *      
!           too high, since the quality of this procedure        *      
!           strongly depends on the machine constant. E.g. on a  *      
!           CONTROL DATA CYBER 175 with a DOUBLE PRECISION       *      
!           machine constant of 2.5D-29, chosing NMAX = 12 has   *      
!           proved adequate at the Computer center of the RWTH   *      
!           Aachen. For other machines a larger MMAX may well be *      
!           meaningful.                                          *      
! LUNIN   : number of the file that contains the nodes for all   *      
!           IRKMs up to order MMAX.                              *      
!           This file can be produced by the SUBROUTINE IRKCOE.  *      
! LUNOUT  : if output of intermediate points is desired (as de-  *      
!           termined by the algorithm), a file with number LUNOUT*      
!           is created for any LUNOUT > 0. If LUNOUT=0, there is *      
!           no output.                                           *      
! LUNPR   : >0; create a log file under number LUNPR, which      *      
!               contains intermediate results and the algorithm  *      
!               flow log.                                        *      
!           =0; no output                                        *      
! EPSM    : DOUBLE PRECISION machine constant                    *      
! EPS     : desired relative accuracy                            *      
! G       : vector G(1:N); the weights G(I) allow varied weighing*      
!                          of the components Y(I) with respect to*      
!                          the machine constant EPS. If all com- *      
!                          ponents are to have the same weight,  *      
!                          then G(I)=1 is chosen for all I.      *      
!                          W A R N I N G:                        *      
!                          If G(I)=0 for one I, division by zero *      
!                          could occur if the corresponding      *      
!                          components of the partial derivatives *      
!                          on the right hand side are equal to   *      
!                          zero. The program will not cover for  *      
!                          this terminal error!                  *      
! X0      : lower limit of the integration interval              *      
! XEND    : upper limit of the integration interval              *      
! Y0      : vector Y0(1:N); initial values Y(X0) at X0           *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS: only used to provide storage space.          *      
! =================                                              *      
! FAC     : vector FAC(1:MMAX); the factoriels FAC(I)=(2*I)!     *      
! DELTAK  : vector DELTAK(0:2*MMAX-1); vector used to determine  *      
!           the theoretically estimated iteration error DK for   *      
!           each iteration step                                  *      
! HEPS    : vector HEPS(1:MMAX); HEPS(M) contains the theoretical*      
!           step size for the IRKM of order M and EPS            *      
! ALPHA   : vector ALPHA(1:MMAX-1);   ) coefficients of          *      
! BETA    : 2-dimensional array       ) the IRKM of              *      
!           BETA(1:MMAX-1,1:MMAX-1);  ) order M                  *      
! A       : vector A(1:MMAX-1);       )                          *      
!                                                                *      
! ALPHAQ  : vector ALPHAQ(1:MMAX);    ) coefficients of          *      
! BETAQ   : 2-dimensional array       ) the IRKM of              *      
!           BETAQ(1:MMAX,1:MMAX);     ) order M + 1              *      
! AQ      : vector AQ(1:MMAX);        )                          *      
!                                                                *      
! DBLEK   : 2-dimensional array DBLEK(1:N,1:MMAX-1); containing  *      
!           the KI of the IRKM of order M                        *      
! DBLEKQ  : 2-dimensional array DBLEKQ(1:N,1:MMAX); contains the *      
!           KI of the IRKM of order M+1                          *      
! DB      : 2-dimensional array DB(1:N,1:MMAX); auxiliary array  *      
!           used to compute the KI for both order IRKMs.         *      
!           this intermediate storage is needed in order to      *      
!           perform the  iteration.                              *      
! Y       : vector Y(1:N); the approximate solution for the      *      
!           method with order M                                  *      
! YOLD    : vector YOLD(1:N); storage of the previous approximate*      
!           solution YQ for accuracy estimate                    *      
! DFDX    : vector DFDX(1:N); derivative DF/DX of the right-hand *      
!           side                                                 *      
! DFDY    : 2-dimensional array DFDY(1:N,1:N); derivative DF/DY  *      
!           of the right-hand side                               *      
! F0      : vector F0(1:N); before each step we evaluate the     *      
!           right-hand side at (X0,Y0)                           *      
! F1      : vector F1(1:N); evaluation of the right-hand side at *      
!           intermediate points, also used when estimating       *      
!           derivatives                                          *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! IERR    : error parameter                                      *      
!           =0; run was successful                               *      
!           =1; no convergence. Possible remedy: increase        *      
!               maximum order MMAX.                              *      
!           =3; too many calls of DES (IFU > MXCALL)             *      
! EPS     : estimate of the largest local relative error         *      
! YQ      : vector YQ(1:N); approximation for the solution of    *      
!           IVP. Obtained by the method of order M+1             *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I, J, K : loop variables                                       *      
! L       : loop variable for the iteration                      *      
! LH      : auxiliary variable for the output of L               *      
! FACLP1  : (L+1)!                                               *      
! IT2     : auxiliary variable I*2                               *      
! M       : currently used order                                 *      
! MM2     : auxiliary variable M*2                               *      
! MM2M1   : auxiliary variable M*2-1                             *      
! MM2M2   : auxiliary variable M*2-2                             *      
! MP1     : auxiliary variable M+1                               *      
! MM2P1   : auxiliary variable M*2+1                             *      
! MMAXM1  : auxiliary variable MMAX-1                            *      
! IFU     : counter for the number of functional evaluations. The*      
!           evaluation of a vector-valued function (N>1) in only *      
!           counted as one evaluation.                           *      
! IPOS    : read counter in the node file, or the most recently  *      
!           read order                                           *      
! IBSP    : ) positioning the read position in the node file     *      
! IGET    : )                                                    *      
! ISC     : counter for the number of intermediate steps         *      
! DELTAG  : difference between two successive approximations     *      
!           YQ(OLD) and YQ(NEW).                                 *      
! DGREL   : relative difference between two successive           *      
!           approximations YQ(OLD) and YQ(NEW).                  *      
! E       : difference between the two approximations Y and YQ   *      
!           for the methods of differing orders.                 *      
! EREL    : relative difference between the two approximations Y *      
!           and YQ for the method with differing orders.         *      
! YQNORM  : L2-norm of YQ                                        *      
! DK      : at first an auxiliry variable for determining        *      
!           DELTAK(L), then used for estimating the iteration    *      
!           error, derived from theoretical considerations in the*      
!           L-th iteration as                                    *      
      SUBROUTINE IMRUKU (DES, N, MMAX, LUNIN, LUNOUT, LUNPR, IERR, EPSM,&
      EPS, G, X0, XEND, Y0, YQ, FAC, DELTAK, HEPS, ALPHA, ALPHAQ, A, AQ,&
      F0, F1, DFDX, Y, YOLD, DFDY, DBLEK, BETA, BETAQ, DBLEKQ, DB)      
!           H**(L+2)/(L+2)! *                                    *      
!              * SQRT (G*(DF/DX*(DF/DY)**L)**2/WGTSUM) * 5,      *      
!           where 5. acts as a safety factor.                    *      
! ICAUSE  : documents the cause for a decrease in step size:     *      
!           =0, no decrease of step size                         *      
!           =1, EREL >= EPS                                      *      
!           =2, DELTAG >= DK                                     *      
!           =3, DGREL >= EPS                                     *      
! MOLD    : order of the previous step. If the new order is equal*      
!           to the old one, the coefficients do not have to be   *      
!           read in again.                                       *      
! IC      : IC serves for dynamic modification of the safety     *      
!           factor. In IC we count the number of step size       *      
!           decreases caused by E<EPS.                           *      
! DELTA   : square root of the machine constant for estimating   *      
!           derivatives with forward difference quotients.       *      
! WGTSUM  : G(1) + G(2) + ... + G(N)                             *      
! SG      : sign for the step size H, this determines the        *      
!           direction of integration.                            *      
! YJ      : auxiliary variable for estimating the partial        *      
!           derivatives DF/DY(J)                                 *      
! YDIFF   : auxiliary variable for determining the difference    *      
! HBEG    : step size at the last step                           *      
! SF      : safety factor, which is used to decrease the step    *      
!           size, it follows from theoretical considerations.    *      
!           SF is modified according to the behaviour of the     *      
!           differential equation while integrating.             *      
! X       : upper bound of the sub-interval of integration from  *      
!           X0 to X0+H=X.                                        *      
! X1      : nodes X0 + H*ALPHA(K) in the interval [X0, X].       *      
! EPSERR  : for safety reasons, we use EPSERR = 10*EPSM as the   *      
!           machine constant on occasion.                        *      
! ERROR   : estimate of the accuracy of the approximate solution *      
! EPSLOC  : estimate of the largest local error                  *      
! SUM     : auxiliary variable for determining the derivatives   *      
!           DF/DX                                                *      
! AW1     : ) variables for determining the amount of work       *      
! AW2     : ) AW(EPS,M) = (N+1+4*M*M)/H(EPS,M)                   *      
! TPOW    : determining of powers of ten                         *      
! ZERO    : logic variable                                       *      
!           =.TRUE., if DF/DX=0 or DF/DY=0                       *      
!           =.FALSE. otherwise                                   *      
! ORDCH   : logic variable for changing the order                *      
!           =.TRUE., if the order has to be increased            *      
!           =.FALSE. otherwise                                   *      
! STOP    : logic variable for stopping the integration          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: IRKRT                                   *      
!                                                                *      
!                                                                *      
!  sources : 1. G. Engeln-Muellges, F. Reutter:                  *      
!               Numerische Mathematik fr Ingenieure, 4th ed.    *      
!               1985.                                            *      
!            2. D. Sommer, see [SOMM67].                         *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 08.29.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!     symbolic constants                                                
!                                                                       
      INTEGER MXCALL, IFU 
!                                                                       
!     Attention: The value 200000 may be too big for an INTEGER-constant
!                                                                       
      PARAMETER (MXCALL = 200000) 
!                                                                       
!     parameters                                                        
!                                                                       
      INTEGER N, MMAX, LUNIN, LUNOUT, LUNPR, IERR 
      DIMENSION G (N), Y0 (N), YQ (N) 
      DOUBLEPRECISION EPSM, EPS, X0, XEND, IRKRT 
!                                                                       
!     parameters (providing storage space)                              
!                                                                       
      DIMENSION FAC (MMAX), F0 (N), F1 (N), DFDX (N), DFDY (N, N),      &
      DELTAK (0:2 * MMAX - 1), HEPS (MMAX), ALPHA (MMAX - 1), BETA (    &
      MMAX - 1, MMAX - 1), A (MMAX - 1), ALPHAQ (MMAX), BETAQ (MMAX,    &
      MMAX), AQ (MMAX), DBLEK (N, MMAX - 1), DBLEKQ (N, MMAX), Y (N),   &
      YOLD (N), DB (N, MMAX)                                            
!                                                                       
!     local variables                                                   
!                                                                       
      INTEGER IT2, I, J, K, L, MM2M1, MM2M2, MP1, IPOS, IBSP, MM2P1,    &
      MM2, ISC, ICAUSE, MOLD, MMAXM1, IC, IGET, LH, M                   
      DOUBLEPRECISION SF, DELTA, WGTSUM, SG, YJ, TPOW, DK, DELTAG, HBEG,&
      X, X1, EPSERR, ERROR, SUM, E, FACLP1, H, YDIFF, EPSLOC, EREL,     &
      DGREL, YQNORM                                                     
      REAL AW1, AW2 
      LOGICAL ZERO, ORDCH, STOP 
!                                                                       
      DATA SF / 0.9D0 / 
!                                                                       
!*****************************************************************      
!*               i n i t i a l i z a t i o n                     *      
!*****************************************************************      
!                                                                       
      IERR = 0 
      EPSERR = 10.0D0 * EPSM 
      EPSLOC = 0.0D0 
      DELTA = DSQRT (EPSM) 
      MMAXM1 = MMAX - 1 
      X = XEND-X0 
      SG = DSIGN (1.0D0, X) 
      IFU = 0 
      ISC = 0 
      M = 0 
      IPOS = 0 
!                                                                       
      REWIND (LUNIN) 
      IF (LUNOUT.GT.0) THEN 
!                                                                       
!        header for output file                                         
!                                                                       
         I = 1 
         WRITE (LUNOUT, 5100) X0, I, Y0 (1) 
         WRITE (LUNOUT, 5200) (I, Y0 (I), I = 2, N) 
         WRITE (LUNOUT, 5300) XEND, EPS, MMAX 
         WRITE (LUNOUT, 7000) LUNPR 
         WRITE (LUNOUT, 7100) 
      ENDIF 
      IF (LUNPR.GT.0) THEN 
!                                                                       
!        header for log file                                            
!                                                                       
         I = 1 
         WRITE (LUNPR, 5100) X0, I, Y0 (1) 
         WRITE (LUNPR, 5200) (I, Y0 (I), I = 2, N) 
         WRITE (LUNPR, 5300) XEND, EPS, MMAX 
         WRITE (LUNPR, 5350) 
         WRITE (LUNPR, 5400) (I, G (I), I = 1, N) 
         WRITE (LUNPR, 6000) LUNOUT 
         WRITE (LUNPR, 6100) 
      ENDIF 
!                                                                       
!     determine the factorials, which                                   
!     are required   FAC(I) = (2*I)!                                    
!                                                                       
      FAC (1) = 2.0D0 
      DO 10 I = 2, MMAX 
         IT2 = I * 2 
         FAC (I) = FAC (I - 1) * DBLE (IT2 * (IT2 - 1) ) 
   10 END DO 
!                                                                       
!     determine the sum of the weights                                  
!                                                                       
      WGTSUM = 0.0D0 
      DO 20 I = 1, N 
         WGTSUM = WGTSUM + G (I) 
   20 END DO 
!                                                                       
!*****************************************************************      
!*                   i n t e g r a t i o n                       *      
!*****************************************************************      
!                                                                       
!     in the following REPEAT-loop, one integration step is             
!     performed in each loop                                            
!                                                                       
 1000 CONTINUE 
      ISC = ISC + 1 
      MOLD = M 
      STOP = .FALSE. 
!                                                                       
!        check number of function calls                                 
!                                                                       
      IF (IFU.GT.MXCALL) THEN 
!                                                                       
!          if necessary stop and set error number                       
!                                                                       
         STOP = .TRUE. 
         IERR = 3 
      ENDIF 
!                                                                       
!        determine approximations for the partial derivatives           
!        of the right-hand side with respect to Y at                    
!        (X0,Y0) with forward difference quotients                      
!                                                                       
      CALL DES (N, X0, Y0, F0) 
      DO 40 J = 1, N 
         YJ = Y0 (J) 
         Y0 (J) = YJ + DELTA 
         CALL DES (N, X0, Y0, F1) 
         DO 30 I = 1, N 
            DFDY (I, J) = (F1 (I) - F0 (I) ) / DELTA 
   30    END DO 
         Y0 (J) = YJ 
   40 END DO 
!                                                                       
!        determine approximations for the partial derivatives           
!        of the right-hand side with respect to X at                    
!        (X0,Y0) with forward difference quotients                      
!                                                                       
      CALL DES (N, X0 + DELTA, Y0, DFDX) 
      DO 50 I = 1, N 
         DFDX (I) = (DFDX (I) - F0 (I) ) / DELTA 
   50 END DO 
      DO 70 I = 1, N 
         SUM = 0.0D0 
         DO 60 J = 1, N 
            SUM = SUM + DFDY (I, J) * F0 (J) 
   60    END DO 
         DFDX (I) = DFDX (I) + SUM 
   70 END DO 
!                                                                       
      IFU = IFU + N + 2 
!                                                                       
!        check whether all partial derivatives with                     
!        respect to Y are equal to 0                                    
!                                                                       
      ZERO = .TRUE. 
      I = 0 
      J = 0 
   80 CONTINUE 
      I = I + 1 
   90 CONTINUE 
      J = J + 1 
      IF (DFDY (I, J) .NE.0.0D0) ZERO = .FALSE. 
      IF (ZERO.AND.J.LT.N) GOTO 90 
      IF (ZERO.AND.I.LT.N) GOTO 80 
!                                                                       
!        check whether all partial derivatives                          
!        with respect to X are equal to 0                               
!                                                                       
      IF (.NOT.ZERO) THEN 
         ZERO = .TRUE. 
         I = 0 
  100    CONTINUE 
         I = I + 1 
         IF (DFDX (I) .NE.0.0D0) ZERO = .FALSE. 
         IF (ZERO.AND.I.LT.N) GOTO 100 
      ENDIF 
!                                                                       
!*****************************************************************      
!*   d e t e r m i n e  t h e  a m o u n t  o f  w o r k                
!*****************************************************************      
!                                                                       
!        estimate the average iteration error, determine                
!        the step size and the amount of work AW                        
!                                                                       
!        the estimated average iteration error at the L-th              
!        iteration is given as                                          
!        H**(L+2)/(L+2)! * SQRT(G*(DF/DX*(DF/DY)**L)**2/WGTSUM).        
!        Here the square root is determined by function IRKRT and it is 
!        stored in DELTAK(L). The iteration error is given as           
!        DK = H**(L+2) * DELTAK(L) / (L+2)! .                           
!                                                                       
!        The step size for the method of order M is HEPS(M)=            
!        EPS*(2*M)! / SQRT(G*(DF/DX*(DF/DY)**(2*M-2))**2/WGTSUM)        
!        **(1/(2*M))                                                    
!                                                                       
!        The amount of work AW is given as AW(EPS,M) =                  
!        (N+1+4*M*M) / HEPS(M). We use the order for which the          
!        amount of work per functional evaluations needed is minimal.   
!        This is the case for the first time for M,                     
!        for which   AW(EPS,M) < AW(EPS,M+1).                           
!                                                                       
      IF (.NOT.ZERO) THEN 
!                                                                       
!           determine the order, the step size, as well as              
!           estimate  the average iteration error if none of the        
!           derivatives is equal to zero                                
!                                                                       
         M = 1 
         DK = 0.0D0 
         DO 110 I = 1, N 
            DK = DK + G (I) * DFDX (I) * DFDX (I) 
  110    END DO 
         DELTAK (0) = DSQRT (DK / WGTSUM) 
         HEPS (1) = DSQRT (EPS * FAC (1) / DELTAK (0) ) 
         AW2 = REAL (N + 5) / REAL (HEPS (1) ) 
         M = 2 
 2000    CONTINUE 
         AW1 = AW2 
         MM2M2 = M * 2 - 2 
         DELTAK (2 * M - 3) = IRKRT (N, DFDY, DFDX, WGTSUM, G, F1) 
         DELTAK (MM2M2) = IRKRT (N, DFDY, DFDX, WGTSUM, G, F1) 
         IF (DELTAK (MM2M2) .NE.0.0D0) THEN 
            HEPS (M) = (EPS * FAC (M) / DELTAK (MM2M2) ) ** (1.0D0 /    &
            DBLE (2 * M) )                                              
            AW2 = REAL (N + 1 + 4 * M * M) / REAL (HEPS (M) ) 
            M = M + 1 
         ELSE 
            ZERO = .TRUE. 
         ENDIF 
         IF (M.LE.MMAX.AND.AW2.LT.AW1.AND..NOT.ZERO) GOTO 2000 
         M = M - 2 
      ENDIF 
      IF (ZERO) THEN 
!                                                                       
!           when determining the step size HEPS(M) we encountered       
!           division by zero. Then we chose the order 3 is chosen and   
!           set DELTAK(L) = EPS * 10**(5-L), for L=1, ..., 5.           
!           For the initial step size  we chose 0.1 .                   
!                                                                       
         M = 3 
         TPOW = 1.0D0 
         DO 120 I = 1, 6 
            DELTAK (6 - I) = EPS * TPOW 
            TPOW = TPOW * 10.0D0 
  120    END DO 
         H = SG * 0.1D0 
      ELSE 
!                                                                       
!           determine the ultimate step size                            
!                                                                       
         H = SG * SF * HEPS (M) 
      ENDIF 
!                                                                       
!        initializations for the new step                               
!                                                                       
      IC = 0 
      X = X0 + H 
      IF (SG * (XEND-X) .LT.0.0D0) THEN 
         H = XEND-X0 
         X = XEND 
      ENDIF 
      ORDCH = .FALSE. 
      IF (M.NE.MOLD) THEN 
!                                                                       
!*****************************************************************      
!*              r e a d   c o e f f i c i e n t s                *      
!*****************************************************************      
!                                                                       
!           if the new order is different from the one in the old       
!           step, the coefficients have to be read in again             
!                                                                       
         IBSP = IPOS - M 
         IGET = M - IPOS 
         IF (IBSP.EQ.0) THEN 
!                                                                       
!              new order is 1 larger than the old one                   
!                                                                       
            BACKSPACE (LUNIN) 
         ELSEIF (IBSP.GT.0) THEN 
!                                                                       
!              new order is smaller than the old one                    
!                                                                       
            DO 130 I = 1, IBSP 
               BACKSPACE (LUNIN) 
  130       END DO 
            BACKSPACE (LUNIN) 
            IGET = 0 
         ELSE 
!                                                                       
!              new order is larger by at least 2                        
!                                                                       
            IGET = IGET - 1 
         ENDIF 
         DO 140 I = 1, IGET 
            READ (LUNIN) 
  140    END DO 
         READ (LUNIN) IPOS, (ALPHA (J), J = 1, IPOS), ( (BETA (J, K),   &
         J = 1, IPOS), K = 1, IPOS), (A (J), J = 1, IPOS)               
         READ (LUNIN) IPOS, (ALPHAQ (J), J = 1, IPOS), ( (BETAQ (J, K), &
         J = 1, IPOS), K = 1, IPOS), (AQ (J), J = 1, IPOS)              
!                                                                       
!           determining frequently used constants                       
!                                                                       
         MP1 = M + 1 
         MM2 = M * 2 
         MM2P1 = MM2 + 1 
         MM2M1 = MM2 - 1 
         MM2M2 = MM2 - 2 
      ENDIF 
!                                                                       
!*****************************************************************      
!*            P E R F O R M I N G  O N E  S T E P                *      
!*****************************************************************      
!                                                                       
!        the following REPEAT-loop is executed until the required       
!        accuracy is achieved, i.e., until the step was executed        
!        successfully. Exception: step width H falls below EPSERR,      
!        which indicates that the procedure does not converge           
!                                                                       
 3000 CONTINUE 
      IF (ORDCH) THEN 
!                                                                       
!              the order is increased by 1. ALPHA, BETA and A           
!              are overwritten with the coefficients ALPHAQ,            
!              BETAQ and AQ                                             
!                                                                       
         ORDCH = .FALSE. 
         M = M + 1 
         MP1 = M + 1 
         DO 160 I = 1, M 
            ALPHA (I) = ALPHAQ (I) 
            DO 150 J = 1, M 
               BETA (J, I) = BETAQ (J, I) 
  150       END DO 
            A (I) = AQ (I) 
  160    END DO 
!                                                                       
!              the new coefficients for ALPHAQ, BETAQ and AQ            
!              are read in                                              
!                                                                       
         READ (LUNIN) IPOS, (ALPHAQ (J), J = 1, IPOS), ( (BETAQ (J, K), &
         J = 1, IPOS), K = 1, IPOS), (AQ (J), J = 1, IPOS)              
!                                                                       
!              determining frequently used constants                    
!                                                                       
         MM2 = M * 2 
         MM2P1 = MM2 + 1 
         MM2M1 = MM2 - 1 
         MM2M2 = MM2 - 2 
         IF (.NOT.ZERO) THEN 
!                                                                       
!                 DELTAK and HEPS have to be determined for the         
!                 new order. The old values may be used in part.        
!                                                                       
            DELTAK (MM2M1) = IRKRT (N, DFDY, DFDX, WGTSUM, G, F1) 
            DELTAK (MM2) = IRKRT (N, DFDY, DFDX, WGTSUM, G, F1) 
            IF (DELTAK (MM2M2) .NE.0.0D0) THEN 
               HEPS (M) = (EPS * FAC (M) / DELTAK (MM2M2) ) ** (1.0D0 / &
               DBLE (MM2) )                                             
               H = SG * SF * HEPS (M) 
            ELSE 
               ZERO = .TRUE. 
            ENDIF 
         ENDIF 
         IF (ZERO) THEN 
!                                                                       
!                 When determining the step size HEPS(M) we have        
!                 encountered division by zero. In this case we set     
!                 DELTAK(L) = EPS * 10**(5-L), for L = 1, ..., 2*M-1    
!                 and chose 0.1 as the initial step size.               
!                                                                       
            TPOW = 1.0D0 
            DO 170 I = 1, MM2 
               DELTAK (MM2 - I) = EPS * TPOW 
               TPOW = TPOW * 10.0D0 
  170       END DO 
            H = SG * 0.1D0 
         ENDIF 
      ENDIF 
!                                                                       
!*****************************************************************      
!*                      i t e r a t i o n                        *      
!*****************************************************************      
!                                                                       
!           Initializing the weights KI and KIQ, i.e., DBLEK and        
!           DBLEKQ to be equal to the value H*F0 for the iteration      
!                                                                       
      DO 190 I = 1, N 
         DO 180 J = 1, M 
            DBLEK (I, J) = H * F0 (I) 
            DBLEKQ (I, J) = H * F0 (I) 
  180    END DO 
         DBLEKQ (I, MP1) = H * F0 (I) 
  190 END DO 
!                                                                       
!           store the old approximation,                                
!           in order to be able to find an error                        
!           estimate after the first iteration                          
!                                                                       
      DO 200 I = 1, N 
         YOLD (I) = Y0 (I) + H * F0 (I) 
  200 END DO 
!                                                                       
!           In the following REPEAT-loop we perform a complete          
!           iteration until the required accuracy is achieved           
!           or until the number of iterations exceeds 2*M+1.            
!           If the latter is the case, the order is increased           
!           and the step is repeated. If a convergence of this          
!           iteration is not occuring, or if the order cannot           
!           be increased further, the whole step is repeated            
!           with a lower step size                                      
!                                                                       
      L = 0 
      FACLP1 = 1.0D0 
 4000 CONTINUE 
!                                                                       
!              check number of function calls                           
!                                                                       
      IF (IFU.GT.MXCALL) THEN 
!                                                                       
!                 if necessary stop and set error number                
!                                                                       
         STOP = .TRUE. 
         IERR = 3 
      ENDIF 
      L = L + 1 
      FACLP1 = FACLP1 * (L + 1) 
      ICAUSE = 0 
!                                                                       
!*****************************************************************      
!*           L-th full iteration to find the weights DBLEK     *        
!*****************************************************************      
!                                                                       
      DO 240 K = 1, M 
!                                                                       
!                 determine intermediate points in the                  
!                 interval [X0, X0+H]                                   
!                                                                       
         X1 = X0 + H * ALPHA (K) 
         DO 220 I = 1, N 
            Y (I) = Y0 (I) 
            DO 210 J = 1, M 
               Y (I) = Y (I) + DBLEK (I, J) * BETA (J, K) 
  210       END DO 
  220    END DO 
!                                                                       
!                 insert the intermediate points into the               
!                 right-hand side of the differential equation          
!                                                                       
         CALL DES (N, X1, Y, F1) 
!                                                                       
!                 perform the full step                                 
!                                                                       
         DO 230 I = 1, N 
            DB (I, K) = H * F1 (I) 
  230    END DO 
  240 END DO 
!                                                                       
!              store the new weights DBLEK after the                    
!              full step                                                
!                                                                       
      DO 260 I = 1, N 
         DO 250 J = 1, M 
            DBLEK (I, J) = DB (I, J) 
  250    END DO 
  260 END DO 
!                                                                       
!*****************************************************************      
!*           L-th iteration for the weights DBLEKQ using a       *      
!*           complete step                                       *      
!*****************************************************************      
!                                                                       
      DO 300 K = 1, MP1 
!                                                                       
!                 determine the intermediate points in the              
!                 interval [X0, X0+H]                                   
!                                                                       
         X1 = X0 + H * ALPHAQ (K) 
         DO 280 I = 1, N 
            YQ (I) = Y0 (I) 
            DO 270 J = 1, MP1 
               YQ (I) = YQ (I) + DBLEKQ (I, J) * BETAQ (J, K) 
  270       END DO 
  280    END DO 
!                                                                       
!                 insert the intermediate points into the               
!                 right-hand side of the differential equation          
!                                                                       
         CALL DES (N, X1, YQ, F1) 
!                                                                       
!                 perform one full step                                 
!                                                                       
         DO 290 I = 1, N 
            DB (I, K) = H * F1 (I) 
  290    END DO 
  300 END DO 
!                                                                       
!              store the new weights DBLEK after the                    
!              complete step                                            
!                                                                       
      DO 320 I = 1, N 
         DO 310 J = 1, MP1 
            DBLEKQ (I, J) = DB (I, J) 
  310    END DO 
  320 END DO 
!                                                                       
      IFU = IFU + MM2P1 
!                                                                       
!*****************************************************************      
!*           approximations  from  the  L-th  iteration          *      
!*****************************************************************      
!                                                                       
!              determine new approximations Y with                      
!              order M, and YQ with order M+1, in the L-th              
!              iteration using RUNGE-KUTTA                              
!                                                                       
      DO 340 I = 1, N 
         Y (I) = Y0 (I) 
         YQ (I) = Y0 (I) 
         DO 330 J = 1, M 
            Y (I) = Y (I) + A (J) * DBLEK (I, J) 
            YQ (I) = YQ (I) + AQ (J) * DBLEKQ (I, J) 
  330    END DO 
         YQ (I) = YQ (I) + AQ (MP1) * DBLEKQ (I, MP1) 
  340 END DO 
!                                                                       
!              determine the absolute and relative average              
!              difference DELTAG and DGREL between two                  
!              successive iterations, as well as the absolute           
!              or relative average difference E or EREL                 
!              between the two approximative solutions                  
!              obtained using methods of different order                
!                                                                       
      DELTAG = 0.0D0 
      E = 0.0D0 
      YQNORM = 0.0D0 
      DO 350 I = 1, N 
         YQNORM = YQNORM + YQ (I) * YQ (I) 
         YDIFF = YQ (I) - YOLD (I) 
         DELTAG = DELTAG + G (I) * YDIFF * YDIFF 
         YDIFF = Y (I) - YQ (I) 
         E = E+G (I) * YDIFF * YDIFF 
         YOLD (I) = YQ (I) 
  350 END DO 
      DELTAG = DSQRT (DELTAG / WGTSUM) 
      E = DSQRT (E / WGTSUM) 
      IF (YQNORM.GT.0.0D0) THEN 
         YQNORM = DSQRT (YQNORM) 
         DGREL = DELTAG / YQNORM 
         EREL = E / YQNORM 
      ELSE 
         DGREL = DELTAG 
         EREL = E 
      ENDIF 
      ERROR = DMAX1 (EREL, DGREL) 
      EPSLOC = DMAX1 (EPSLOC, ERROR) 
!                                                                       
!*****************************************************************      
!*             test  for  stopping the  iteration                *      
!*****************************************************************      
!                                                                       
      LH = L 
      IF (EREL.GE.EPS) THEN 
!                                                                       
!                 the difference between two approximations for         
!                 different orders differ by a term of the same         
!                 order 2*M+1 when using full steps as that of the      
!                 approximate solution.                                 
!                 If this is not the case the step size is decreased    
!                 according to the theory.                              
!                                                                       
         ICAUSE = 1 
         IF (IC.NE.0) SF = SF * 0.9D0 
         IC = IC + 1 
         HBEG = H 
         H = H * SF * (0.5D0 * EPS / E) ** (1.0D0 / DBLE (MM2P1) ) 
         L = 0 
      ELSE 
!                                                                       
!                 compute the estimated average iteration error         
!                                                                       
         IF (ZERO) THEN 
            DK = DELTAK (L - 1) 
         ELSE 
            DK = 5.0D0 * H** (L + 1) * DELTAK (L - 1) / FACLP1 
         ENDIF 
!                                                                       
!                 When computing H only the local error is considered,  
!                 hence we need to test for convergence of the          
!                 iteration procedure                                   
!                                                                       
         IF (DELTAG.GE.DK) THEN 
            ICAUSE = 2 
            HBEG = H 
            H = 0.6D0 * H 
            SF = 0.8D0 * SF 
            L = 0 
         ELSE 
            IF (DGREL.GE.EPS.AND.M.GE.MMAXM1.AND.L.GE.MM2M1) THEN 
!                                                                       
!                       if the maximum number of iterations is          
!                       reached but without the desired accuracy,       
!                       then the order should be increased. If this     
!                       is not possible, since e.g. no more nodes       
!                       are available, we attempt to come to a successfu
!                       end by decreasing the step size                 
!                                                                       
               ICAUSE = 3 
               HBEG = H 
               H = 0.8D0 * H 
               L = 0 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (L.EQ.0) THEN 
!                                                                       
!                 L=0 indicates that the step  cannot be                
!                 completed successfully                                
!                                                                       
         IF (LUNPR.GT.0) THEN 
            WRITE (LUNPR, 6200) ISC, M, HBEG, X, YQ (1), ERROR, IFU, LH,&
            ICAUSE                                                      
            WRITE (LUNPR, 6300) (YQ (I), I = 2, N) 
         ENDIF 
         IF (H.LT.EPSERR) THEN 
!                                                                       
!                    the procedure does not converge                    
!                                                                       
            STOP = .TRUE. 
            IERR = 1 
         ELSE 
!                                                                       
!                    the calculations of the last step are              
!                    cancelled and the step is repeated                 
!                                                                       
            FACLP1 = 1.0D0 
            X = X0 + H 
            IF (SG * (XEND-X) .LT.0.0D0) THEN 
               H = XEND-X0 
               X = XEND 
            ENDIF 
            DO 370 I = 1, N 
               DO 360 J = 1, M 
                  DBLEK (I, J) = H * F0 (I) 
                  DBLEKQ (I, J) = H * F0 (I) 
  360          END DO 
               DBLEKQ (I, MP1) = H * F0 (I) 
  370       END DO 
            DO 380 I = 1, N 
               YOLD (I) = Y0 (I) + H * F0 (I) 
  380       END DO 
         ENDIF 
      ENDIF 
!                                                                       
!*****************************************************************      
!*                l o o p    i n q u i r i e s                   *      
!*****************************************************************      
!                                                                       
      IF ( ( (L.LT.MM2M1.AND.DGREL.GE.EPS) .OR.L.EQ.0) .AND..NOT.STOP)  &
      GOTO 4000                                                         
      IF (DGREL.GE.EPS.AND.MP1.LT.MMAX.AND..NOT.STOP) THEN 
!                                                                       
!              the desired accuracy has not been reached after          
!              the theoretically maximal number of iterations. The      
!              order is increased and the step is repeated              
!                                                                       
         ORDCH = .TRUE. 
         ICAUSE = 4 
         IF (LUNPR.GT.0) THEN 
            WRITE (LUNPR, 6200) ISC, M, H, X, YQ (1), ERROR, IFU, L,    &
            ICAUSE                                                      
            WRITE (LUNPR, 6300) (YQ (I), I = 2, N) 
         ENDIF 
      ELSE 
!                                                                       
!              prevent from overstep the maximal order                  
!                                                                       
         ORDCH = .FALSE. 
      ENDIF 
      IF (ORDCH) GOTO 3000 
      IF (.NOT.STOP) THEN 
!                                                                       
!           the step was successfully completed                         
!                                                                       
         IF (LUNPR.GT.0) THEN 
            WRITE (LUNPR, 6200) ISC, M, H, X, YQ (1), ERROR, IFU, L,    &
            ICAUSE                                                      
            WRITE (LUNPR, 6300) (YQ (I), I = 2, N) 
         ENDIF 
         IF (LUNOUT.GT.0) THEN 
            I = 1 
            WRITE (LUNOUT, 7200) ISC, X, I, YQ (1), ERROR 
            WRITE (LUNOUT, 7300) (I, YQ (I), I = 2, N) 
         ENDIF 
      ENDIF 
!                                                                       
      IF (X.NE.XEND.AND..NOT.STOP) THEN 
!                                                                       
!           the end of the integration interval has not                 
!           been reached                                                
!                                                                       
         IF (IC.GT.1) SF = SF / 0.97D0 
         X0 = X 
         DO 390 I = 1, N 
            Y0 (I) = YQ (I) 
  390    END DO 
      ELSE 
!                                                                       
!           program stop since the initial value problem was            
!           solved (IERR=0), or because the procedure does not          
!           converge (IERR=1)                                           
!                                                                       
         STOP = .TRUE. 
         EPS = EPSLOC 
      ENDIF 
      IF (.NOT.STOP) GOTO 1000 
!                                                                       
!*****************************************************************      
!*           f o r m a t  s t a t e m e n t s                    *      
!*****************************************************************      
!                                                                       
 5100 FORMAT ('1','INITIAL CONDITION:',//,                              &
     &        1X,10X,'X0',10X,2X,'COMP.',2X,10X,'Y0',/,                 &
     &        1X,E22.15,      1X,I4,1X  ,2X,E22.15)                     
 5200 FORMAT (1X,            23X,I4,1X  ,2X,E22.15) 
 5300 FORMAT ('0','RIGHT ENDPOINT OF THE INTERVAL OF INTEGRATION: ',    &
     &        E22.15,/,30X,'DESIRED ACCURACY: ',E22.15,/                &
     &         ,25X,'LARGEST ORDER ALLOWED:',I4)                        
 5350 FORMAT ('0','REASON FOR DECREASE IN STEP SIZE' ,/,5X              &
     &        ,'0  NO DECREASE',/,5X                                    &
     &        ,'1  EREL >= EPS',/,5X,'2  DELTAG >= DK',/,5X             &
     &        ,'3  DGREL >= EPS')                                       
 5400 FORMAT ('0','WEIGHTS G:       COMP.  ',10X,'G',/,                 &
     &         (17X,I4,3X,E22.15))                                      
 6000 FORMAT ('0','NUMBER OF THE OUTPUT FILE (0=NONE):',I4) 
 6100 FORMAT ('0',' STEP  ORDER STEP SIZE  UPPER BOUND APPROXIMATION',  &
     &        ' ERROR      FUNCTION ITERATION CAUSE',/,18X,'H',24X,     &
     &        'Y',7X,'ESTIMATE',3X,'CALLS',4X,'STEPS')                  
 6200 FORMAT (1X,I6,2X,I3,2X,E10.3,1X,E10.3,2X,E13.6,1X,E10.3,1X,       &
     &        I8,1X,I9,3X,I1)                                           
 6300 FORMAT (36X,E13.6) 
 7000 FORMAT ('0','NUMBER OF THE LOG FILE (0=NONE):',I4) 
 7100 FORMAT ('0',' STEP  BOUND FOR UPPER LIMIT  COMP.',2X,             &
     &        'APPROXIMATE SOLUTION     ERROR ESTIMATE',/,7X,           &
     &        ' OF INTEGRATIONINTERVAL ')                               
 7200 FORMAT (1X,I6,1X,E22.15,1X,I4,2X,D24.17,1X,E18.11) 
 7300 FORMAT (31X,I4,2X,D24.17) 
!                                                                       
      RETURN 
      END SUBROUTINE IMRUKU                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION IRKRT (N, DFDY, DFDX, WGTSUM, G, DFDXH) 
!                                                                       
!*****************************************************************      
!                                                                *      
! The DOUBLE PRECISION FUNCTION IRKRT determines square roots,   *      
! that occur in SUBROUTINE IMRUKU when determining the step size *      
! H and the estimated average iteration error DK or DELTAK(L),   *      
! L=0, 1, ...., 2*M+1.                                           *      
!                                                                *      
! IRKRT = SQRT(G*(DF/DX*(DF/DY)**K))**2/WGTSUM), for K=1,..,2*M-1*      
!                                                                *      
! Here we determine (DF/DY)**K during each call of IRKRT, by     *      
! storing (DF/DY)**(K-1)*DF/DX in vector previously used for     *      
! DF/DX and DFDX.                                                *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS: (compare description of variables in IMRUKU) *      
! =================                                              *      
! N       : dimension                                            *      
! DFDY    : 2-dimensional array DFDY(1:N,1:N); the derivative    *      
!           DF/DY of the right-hand side of the differential     *      
!           equation with respect to Y from IMRUKU at (X0,Y0).   *      
! DFDX    : vector DFDX(1:N); the derivative DF/DX of the right- *      
!           hand side of the differential equation with respect  *      
!           to X from IMRUKU or (DF/DY)**(K-1) * DF/DX at (X0,Y0)*      
! G       : vector G(1:N) of weights                             *      
! WGTSUM  : G(1) + G(2) + ... + G(N)                             *      
! DFDXH   : vector DFDXH(1:N); auxiliary storage                 *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! IRKRT   : computed square root (compare program description)   *      
! DFDX    : vector DFDX(1:N);  (DF/DY)**K * DF/DX                *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I,J     : control variables                                    *      
! SUM     : scalar products DFDY(I,..)*DFDX                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 08.12.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!     parameter                                                         
!                                                                       
      INTEGER N 
      DIMENSION DFDX (N), DFDY (N, N), G (N), DFDXH (N) 
      DOUBLEPRECISION WGTSUM 
!                                                                       
!     local variables                                                   
!                                                                       
      INTEGER I, J 
      DOUBLEPRECISION SUM 
!                                                                       
      DO 20 I = 1, N 
         SUM = 0.0D0 
         DO 10 J = 1, N 
            SUM = SUM + DFDY (I, J) * DFDX (J) 
   10    END DO 
         DFDXH (I) = SUM 
   20 END DO 
      SUM = 0.0D0 
      DO 30 J = 1, N 
         DFDX (J) = DFDXH (J) 
         SUM = SUM + DFDX (J) **2 * G (J) 
   30 END DO 
      IRKRT = DSQRT (SUM / WGTSUM) 
!                                                                       
      RETURN 
      END FUNCTION IRKRT                            

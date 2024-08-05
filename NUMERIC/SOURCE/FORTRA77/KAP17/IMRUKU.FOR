C[BA*)
C[LE*)
      SUBROUTINE IMRUKU (DES,N,MMAX,LUNIN,LUNOUT,LUNPR,IERR,EPSM,
C[IX{IMRUKU}*)
     1                   EPS,G,X0,XEND,Y0,YQ,
     2                   FAC,DELTAK,HEPS,ALPHA,ALPHAQ,A,AQ,
     3                   F0,F1,DFDX,Y,YOLD,DFDY,
     4                   DBLEK,BETA,BETAQ,DBLEKQ,DB)
C
C*****************************************************************
C                                                                *
C IMRUKU solves an initial value problem (IVP) for systems of N  *
C given differential equations of 1st order.                     *
C The solution of the IVP is obtained with an implicit RUNGE-    *
C KUTTA method (IRKM), where we use step size control and control*
C the order of the IRKM as well.                                 *
C[BE*)
C For this the subroutine IMRUKU needs IRKMs of orders 1 to MMAX.*
C The coefficients that define these implicit RUNGE-KUTTA methods*
C are assumed to be available in a file available at the input   *
C LUNIN. This file can be created using the SUBROUTINE IRKCOE.   *
C                                                                *
C Before each RUNGE-KUTTA step we precheck the efficiency AW in  *
C order to find the optimal order in relation to the number of   *
C functional evaluations according to                            *
C         AW(EPS,M) = (N+1+4*M**2)/H(EPS,M).                     *
C Here H(EPS,M) is the step size for the order M and the accuracy*
C bound EPS. For each step we chose an M for which the work AW   *
C required for the IRKMs of order 1 to M+1 satisfy:              *
C        AW(EPS,I) >  AW(EPS,I+1) for I=1, ..., M-1              *
C  and   AW(EPS,M) <= AW(EPS,M+1);                               *
C  or    M=MMAX-1, if there is no suitable M between 1 and MMAX-1*
C                                                                *
C For this M we chose a step size H which usually is taken equal *
C to the theoretical step size H(EPS,M).                         *
C The RUNGE-KUTTA step is performed for two IRKMs. One IRKM has  *
C the order M with coefficient A, ALPHA and BETA. The other IRKM *
C is of order M+1 for the coefficients AQ, ALPHAQ and BETAQ.     *
C Determining the KI, called DBLEK and DBLEKQ in the program, is *
C done iteratively.                                              *
C                                                                *
C If in an iteration step, the relative difference EREL of the   *
C two approximations Y and YQ, that were obtained by two         *
C different methods, exceeds EPS, then the step size H is reduced*
C to H * SF * (0.5 * EPS * E)**(1/(2*M+1)). Here SF is a safety  *
C factor and E is the absolute difference of Y and YQ. Then the  *
C last step is repeated with the new step size H.                *
C If in an iteration step, the difference between two successive *
C approximations YQ(ALT) and YQ(NEU) exceeds the estimated       *
C theoretical iteration error DK, no convergence can be expected.*
C In this case the step size is decreased to 6/10 of the old     *
C step size and the last integration step is repeated.           *
C Iteration is continued until the relative difference between   *
C two successive approximations becomes smaller than the required*
C accuracy EPS.                                                  *
C Based on theoretical considerations this should be the case    *
C after at most 2*M+1 iterations.                                *
C However, if this does not happen, the order is increased to    *
C M + 1 and the last step is recomputed. If this would exceed    *
C MMAX-1 we alter the step size to 0.8*H instead.                *
C If, due to our step size control, H becomes less than ten      *
C times the machine constant, we stop : no convergence can be    *
C expected on the computer used.                                 *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C DES     : name of a subroutine, that represents the system of  *
C           differential equations. This subroutine has to be    *
C           declared as EXTERNAL in the calling program. It has  *
C           the form:                                            *
C                 SUBROUTINE DES (N,X,Y,F)                       *
C                 INTEGER N                                      *
C                 DOUBLE PRECISION X, Y(N), F(N)                 *
C                 F(1) = F1(X,Y(1),Y(2),...,Y(N))                *
C                 F(2) = F2(X,Y(1),Y(2),...,Y(N))                *
C                 -------------------------------                *
C                 F(N) = FN(X,Y(1),Y(2),...,Y(N))                *
C                 RETURN                                         *
C                 END                                            *
C N       : The number of differential equations                 *
C MMAX    : maximally allowed order, i.e. the highest order for  *
C           which coefficients are available in an auxiliary     *
C           file addressed at LUNIN. MMAX has to be >= 5.        *
C           It is not advisable to chose the maximum order NMAX  *
C           too high, since the quality of this procedure        *
C           strongly depends on the machine constant. E.g. on a  *
C           CONTROL DATA CYBER 175 with a DOUBLE PRECISION       *
C           machine constant of 2.5D-29, chosing NMAX = 12 has   *
C           proved adequate at the Computer center of the RWTH   *
C           Aachen. For other machines a larger MMAX may well be *
C           meaningful.                                          *
C LUNIN   : number of the file that contains the nodes for all   *
C           IRKMs up to order MMAX.                              *
C           This file can be produced by the SUBROUTINE IRKCOE.  *
C LUNOUT  : if output of intermediate points is desired (as de-  *
C           termined by the algorithm), a file with number LUNOUT*
C           is created for any LUNOUT > 0. If LUNOUT=0, there is *
C           no output.                                           *
C LUNPR   : >0; create a log file under number LUNPR, which      *
C               contains intermediate results and the algorithm  *
C               flow log.                                        *
C           =0; no output                                        *
C EPSM    : DOUBLE PRECISION machine constant                    *
C EPS     : desired relative accuracy                            *
C G       : vector G(1:N); the weights G(I) allow varied weighing*
C                          of the components Y(I) with respect to*
C                          the machine constant EPS. If all com- *
C                          ponents are to have the same weight,  *
C                          then G(I)=1 is chosen for all I.      *
C                          W A R N I N G:                        *
C                          If G(I)=0 for one I, division by zero *
C                          could occur if the corresponding      *
C                          components of the partial derivatives *
C                          on the right hand side are equal to   *
C                          zero. The program will not cover for  *
C                          this terminal error!                  *
C X0      : lower limit of the integration interval              *
C XEND    : upper limit of the integration interval              *
C Y0      : vector Y0(1:N); initial values Y(X0) at X0           *
C                                                                *
C                                                                *
C INPUT PARAMETERS: only used to provide storage space.          *
C =================                                              *
C FAC     : vector FAC(1:MMAX); the factoriels FAC(I)=(2*I)!     *
C DELTAK  : vector DELTAK(0:2*MMAX-1); vector used to determine  *
C           the theoretically estimated iteration error DK for   *
C           each iteration step                                  *
C HEPS    : vector HEPS(1:MMAX); HEPS(M) contains the theoretical*
C           step size for the IRKM of order M and EPS            *
C ALPHA   : vector ALPHA(1:MMAX-1);   ) coefficients of          *
C BETA    : 2-dimensional array       ) the IRKM of              *
C           BETA(1:MMAX-1,1:MMAX-1);  ) order M                  *
C A       : vector A(1:MMAX-1);       )                          *
C                                                                *
C ALPHAQ  : vector ALPHAQ(1:MMAX);    ) coefficients of          *
C BETAQ   : 2-dimensional array       ) the IRKM of              *
C           BETAQ(1:MMAX,1:MMAX);     ) order M + 1              *
C AQ      : vector AQ(1:MMAX);        )                          *
C                                                                *
C DBLEK   : 2-dimensional array DBLEK(1:N,1:MMAX-1); containing  *
C           the KI of the IRKM of order M                        *
C DBLEKQ  : 2-dimensional array DBLEKQ(1:N,1:MMAX); contains the *
C           KI of the IRKM of order M+1                          *
C DB      : 2-dimensional array DB(1:N,1:MMAX); auxiliary array  *
C           used to compute the KI for both order IRKMs.         *
C           this intermediate storage is needed in order to      *
C           perform the  iteration.                              *
C Y       : vector Y(1:N); the approximate solution for the      *
C           method with order M                                  *
C YOLD    : vector YOLD(1:N); storage of the previous approximate*
C           solution YQ for accuracy estimate                    *
C DFDX    : vector DFDX(1:N); derivative DF/DX of the right-hand *
C           side                                                 *
C DFDY    : 2-dimensional array DFDY(1:N,1:N); derivative DF/DY  *
C           of the right-hand side                               *
C F0      : vector F0(1:N); before each step we evaluate the     *
C           right-hand side at (X0,Y0)                           *
C F1      : vector F1(1:N); evaluation of the right-hand side at *
C           intermediate points, also used when estimating       *
C           derivatives                                          *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C IERR    : error parameter                                      *
C           =0; run was successful                               *
C           =1; no convergence. Possible remedy: increase        *
C               maximum order MMAX.                              *
C           =3; too many calls of DES (IFU > MXCALL)             *
C EPS     : estimate of the largest local relative error         *
C YQ      : vector YQ(1:N); approximation for the solution of    *
C           IVP. Obtained by the method of order M+1             *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C I, J, K : loop variables                                       *
C L       : loop variable for the iteration                      *
C LH      : auxiliary variable for the output of L               *
C FACLP1  : (L+1)!                                               *
C IT2     : auxiliary variable I*2                               *
C M       : currently used order                                 *
C MM2     : auxiliary variable M*2                               *
C MM2M1   : auxiliary variable M*2-1                             *
C MM2M2   : auxiliary variable M*2-2                             *
C MP1     : auxiliary variable M+1                               *
C MM2P1   : auxiliary variable M*2+1                             *
C MMAXM1  : auxiliary variable MMAX-1                            *
C IFU     : counter for the number of functional evaluations. The*
C           evaluation of a vector-valued function (N>1) in only *
C           counted as one evaluation.                           *
C IPOS    : read counter in the node file, or the most recently  *
C           read order                                           *
C IBSP    : ) positioning the read position in the node file     *
C IGET    : )                                                    *
C ISC     : counter for the number of intermediate steps         *
C DELTAG  : difference between two successive approximations     *
C           YQ(OLD) and YQ(NEW).                                 *
C DGREL   : relative difference between two successive           *
C           approximations YQ(OLD) and YQ(NEW).                  *
C E       : difference between the two approximations Y and YQ   *
C           for the methods of differing orders.                 *
C EREL    : relative difference between the two approximations Y *
C           and YQ for the method with differing orders.         *
C YQNORM  : L2-norm of YQ                                        *
C DK      : at first an auxiliry variable for determining        *
C           DELTAK(L), then used for estimating the iteration    *
C           error, derived from theoretical considerations in the*
C           L-th iteration as                                    *
C           H**(L+2)/(L+2)! *                                    *
C              * SQRT (G*(DF/DX*(DF/DY)**L)**2/WGTSUM) * 5,      *
C           where 5. acts as a safety factor.                    *
C ICAUSE  : documents the cause for a decrease in step size:     *
C           =0, no decrease of step size                         *
C           =1, EREL >= EPS                                      *
C           =2, DELTAG >= DK                                     *
C           =3, DGREL >= EPS                                     *
C MOLD    : order of the previous step. If the new order is equal*
C           to the old one, the coefficients do not have to be   *
C           read in again.                                       *
C IC      : IC serves for dynamic modification of the safety     *
C           factor. In IC we count the number of step size       *
C           decreases caused by E<EPS.                           *
C DELTA   : square root of the machine constant for estimating   *
C           derivatives with forward difference quotients.       *
C WGTSUM  : G(1) + G(2) + ... + G(N)                             *
C SG      : sign for the step size H, this determines the        *
C           direction of integration.                            *
C YJ      : auxiliary variable for estimating the partial        *
C           derivatives DF/DY(J)                                 *
C YDIFF   : auxiliary variable for determining the difference    *
C HBEG    : step size at the last step                           *
C SF      : safety factor, which is used to decrease the step    *
C           size, it follows from theoretical considerations.    *
C           SF is modified according to the behaviour of the     *
C           differential equation while integrating.             *
C X       : upper bound of the sub-interval of integration from  *
C           X0 to X0+H=X.                                        *
C X1      : nodes X0 + H*ALPHA(K) in the interval [X0, X].       *
C EPSERR  : for safety reasons, we use EPSERR = 10*EPSM as the   *
C           machine constant on occasion.                        *
C ERROR   : estimate of the accuracy of the approximate solution *
C EPSLOC  : estimate of the largest local error                  *
C SUM     : auxiliary variable for determining the derivatives   *
C           DF/DX                                                *
C AW1     : ) variables for determining the amount of work       *
C AW2     : ) AW(EPS,M) = (N+1+4*M*M)/H(EPS,M)                   *
C TPOW    : determining of powers of ten                         *
C ZERO    : logic variable                                       *
C           =.TRUE., if DF/DX=0 or DF/DY=0                       *
C           =.FALSE. otherwise                                   *
C ORDCH   : logic variable for changing the order                *
C           =.TRUE., if the order has to be increased            *
C           =.FALSE. otherwise                                   *
C STOP    : logic variable for stopping the integration          *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: IRKRT                                   *
C                                                                *
C                                                                *
C  sources : 1. G. Engeln-Muellges, F. Reutter:                  *
C               Numerische Mathematik fr Ingenieure, 4th ed.    *
C               1985.                                            *
C            2. D. Sommer, see [SOMM67].                         *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Thomas Eul                                         *
C  date     : 08.29.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     symbolic constants
C
      INTEGER*4 MXCALL,IFU
C
C     Attention: The value 200000 may be too big for an INTEGER-constant
C
      PARAMETER (MXCALL = 200000)
C
C     parameters
C
      INTEGER   N, MMAX, LUNIN, LUNOUT, LUNPR, IERR
      DIMENSION G(N), Y0(N), YQ(N)
      DOUBLE PRECISION EPSM, EPS, X0, XEND, IRKRT
C
C     parameters (providing storage space)
C
      DIMENSION FAC(MMAX), F0(N), F1(N), DFDX(N), DFDY(N,N),
     1          DELTAK(0:2*MMAX-1), HEPS(MMAX),
     2          ALPHA(MMAX-1), BETA(MMAX-1,MMAX-1), A(MMAX-1),
     3          ALPHAQ(MMAX),  BETAQ(MMAX,MMAX),    AQ(MMAX),
     4          DBLEK(N,MMAX-1), DBLEKQ(N,MMAX),
     5          Y(N), YOLD(N), DB(N,MMAX)
C
C     local variables
C
      INTEGER IT2, I, J, K, L, MM2M1, MM2M2, MP1, IPOS, IBSP,
     1        MM2P1, MM2, ISC, ICAUSE, MOLD, MMAXM1, IC, IGET, LH,
     2        M
      DOUBLE PRECISION SF, DELTA, WGTSUM, SG, YJ, TPOW, DK, DELTAG
     1                ,HBEG, X, X1, EPSERR, ERROR, SUM, E, FACLP1
     2                ,H, YDIFF, EPSLOC, EREL, DGREL, YQNORM
      REAL AW1, AW2
      LOGICAL ZERO, ORDCH, STOP
C
      DATA SF/0.9D0/
C
C*****************************************************************
C*               i n i t i a l i z a t i o n                     *
C*****************************************************************
C
      IERR   = 0
      EPSERR = 10.0D0 * EPSM
      EPSLOC = 0.0D0
      DELTA  = DSQRT(EPSM)
      MMAXM1 = MMAX - 1
      X      = XEND - X0
      SG     = DSIGN(1.0D0,X)
      IFU    = 0
      ISC    = 0
      M      = 0
      IPOS   = 0
C
      REWIND (LUNIN)
      IF (LUNOUT .GT. 0) THEN
C
C        header for output file
C
         I = 1
         WRITE (LUNOUT,5100) X0,I,Y0(1)
         WRITE (LUNOUT,5200) (I,Y0(I),I=2,N)
         WRITE (LUNOUT,5300) XEND,EPS,MMAX
         WRITE (LUNOUT,7000) LUNPR
         WRITE (LUNOUT,7100)
      ENDIF
      IF (LUNPR .GT. 0) THEN
C
C        header for log file
C
         I = 1
         WRITE (LUNPR,5100) X0,I,Y0(1)
         WRITE (LUNPR,5200) (I,Y0(I),I=2,N)
         WRITE (LUNPR,5300) XEND,EPS,MMAX
         WRITE (LUNPR,5350)
         WRITE (LUNPR,5400) (I,G(I),I=1,N)
         WRITE (LUNPR,6000) LUNOUT
         WRITE (LUNPR,6100)
      ENDIF
C
C     determine the factorials, which
C     are required   FAC(I) = (2*I)!
C
      FAC(1) = 2.0D0
      DO 10 I = 2,MMAX
         IT2 = I * 2
         FAC(I) = FAC(I-1) * DBLE(IT2*(IT2-1))
   10 CONTINUE
C
C     determine the sum of the weights
C
      WGTSUM = 0.0D0
      DO 20 I = 1,N
         WGTSUM = WGTSUM + G(I)
   20 CONTINUE
C
C*****************************************************************
C*                   i n t e g r a t i o n                       *
C*****************************************************************
C
C     in the following REPEAT-loop, one integration step is
C     performed in each loop
C
 1000 CONTINUE
         ISC  = ISC + 1
         MOLD = M
         STOP = .FALSE.
C
C        check number of function calls
C
         IF (IFU .GT. MXCALL) THEN
C
C          if necessary stop and set error number
C
            STOP = .TRUE.
            IERR = 3
         ENDIF
C
C        determine approximations for the partial derivatives
C        of the right-hand side with respect to Y at
C        (X0,Y0) with forward difference quotients
C
         CALL DES (N,X0,Y0,F0)
         DO 40 J = 1,N
            YJ    = Y0(J)
            Y0(J) = YJ + DELTA
            CALL DES (N,X0,Y0,F1)
            DO 30 I = 1,N
               DFDY(I,J) = (F1(I)-F0(I))/DELTA
   30       CONTINUE
            Y0(J) = YJ
   40    CONTINUE
C
C        determine approximations for the partial derivatives
C        of the right-hand side with respect to X at
C        (X0,Y0) with forward difference quotients
C
         CALL DES (N,X0+DELTA,Y0,DFDX)
         DO 50 I = 1,N
            DFDX(I) = (DFDX(I)-F0(I))/DELTA
   50    CONTINUE
         DO 70 I = 1,N
            SUM = 0.0D0
            DO 60 J = 1,N
               SUM = SUM + DFDY(I,J)*F0(J)
   60       CONTINUE
            DFDX(I) = DFDX(I) + SUM
   70    CONTINUE
C
         IFU = IFU + N + 2
C
C        check whether all partial derivatives with
C        respect to Y are equal to 0
C
         ZERO = .TRUE.
         I    = 0
         J    = 0
   80    CONTINUE
            I = I + 1
   90       CONTINUE
               J = J + 1
               IF (DFDY(I,J) .NE. 0.0D0) ZERO = .FALSE.
            IF (ZERO .AND. J .LT. N) GO TO 90
         IF (ZERO .AND. I .LT. N) GO TO 80
C
C        check whether all partial derivatives
C        with respect to X are equal to 0
C
         IF (.NOT. ZERO) THEN
            ZERO = .TRUE.
            I = 0
  100       CONTINUE
               I = I + 1
               IF (DFDX(I) .NE. 0.0D0) ZERO = .FALSE.
            IF (ZERO .AND. I .LT. N) GO TO 100
         ENDIF
C
C*****************************************************************
C*   d e t e r m i n e  t h e  a m o u n t  o f  w o r k
C*****************************************************************
C
C        estimate the average iteration error, determine
C        the step size and the amount of work AW
C
C        the estimated average iteration error at the L-th
C        iteration is given as
C        H**(L+2)/(L+2)! * SQRT(G*(DF/DX*(DF/DY)**L)**2/WGTSUM).
C        Here the square root is determined by function IRKRT and it is
C        stored in DELTAK(L). The iteration error is given as
C        DK = H**(L+2) * DELTAK(L) / (L+2)! .
C
C        The step size for the method of order M is HEPS(M)=
C        EPS*(2*M)! / SQRT(G*(DF/DX*(DF/DY)**(2*M-2))**2/WGTSUM)
C        **(1/(2*M))
C
C        The amount of work AW is given as AW(EPS,M) =
C        (N+1+4*M*M) / HEPS(M). We use the order for which the
C        amount of work per functional evaluations needed is minimal.
C        This is the case for the first time for M,
C        for which   AW(EPS,M) < AW(EPS,M+1).
C
         IF (.NOT. ZERO) THEN
C
C           determine the order, the step size, as well as
C           estimate  the average iteration error if none of the
C           derivatives is equal to zero
C
            M  = 1
            DK = 0.0D0
            DO 110 I = 1,N
               DK = DK + G(I)*DFDX(I)*DFDX(I)
  110       CONTINUE
            DELTAK(0) = DSQRT(DK/WGTSUM)
            HEPS(1)   = DSQRT (EPS * FAC(1) / DELTAK(0))
            AW2       = REAL(N+5) / REAL(HEPS(1))
            M         = 2
 2000       CONTINUE
               AW1           = AW2
               MM2M2         = M*2 - 2
               DELTAK(2*M-3) = IRKRT (N,DFDY,DFDX,WGTSUM,G,F1)
               DELTAK(MM2M2) = IRKRT (N,DFDY,DFDX,WGTSUM,G,F1)
               IF (DELTAK(MM2M2) .NE. 0.0D0) THEN
                  HEPS(M) = (EPS * FAC(M) / DELTAK(MM2M2))
     1                       ** (1.0D0 / DBLE(2*M))
                  AW2     = REAL(N+1+4*M*M) / REAL(HEPS(M))
                  M       = M + 1
               ELSE
                  ZERO = .TRUE.
               ENDIF
            IF (M .LE. MMAX .AND. AW2 .LT. AW1 .AND.
     1                                         .NOT. ZERO) GOTO 2000
            M = M - 2
         ENDIF
         IF (ZERO) THEN
C
C           when determining the step size HEPS(M) we encountered
C           division by zero. Then we chose the order 3 is chosen and
C           set DELTAK(L) = EPS * 10**(5-L), for L=1, ..., 5.
C           For the initial step size  we chose 0.1 .
C
            M    = 3
            TPOW = 1.0D0
            DO 120 I = 1,6
               DELTAK(6-I) = EPS * TPOW
               TPOW        = TPOW * 10.0D0
  120       CONTINUE
            H = SG * 0.1D0
         ELSE
C
C           determine the ultimate step size
C
            H = SG * SF * HEPS(M)
         ENDIF
C
C        initializations for the new step
C
         IC = 0
         X  = X0 + H
         IF (SG*(XEND-X) .LT. 0.0D0) THEN
            H = XEND - X0
            X = XEND
         ENDIF
         ORDCH = .FALSE.
         IF (M .NE. MOLD) THEN
C
C*****************************************************************
C*              r e a d   c o e f f i c i e n t s                *
C*****************************************************************
C
C           if the new order is different from the one in the old
C           step, the coefficients have to be read in again
C
            IBSP = IPOS - M
            IGET = M - IPOS
            IF (IBSP .EQ. 0) THEN
C
C              new order is 1 larger than the old one
C
               BACKSPACE (LUNIN)
            ELSEIF (IBSP .GT. 0) THEN
C
C              new order is smaller than the old one
C
               DO 130 I = 1,IBSP
                  BACKSPACE (LUNIN)
  130          CONTINUE
               BACKSPACE (LUNIN)
               IGET = 0
            ELSE
C
C              new order is larger by at least 2
C
               IGET = IGET - 1
            ENDIF
            DO 140 I = 1,IGET
               READ (LUNIN)
  140       CONTINUE
            READ (LUNIN) IPOS,(ALPHA(J),J=1,IPOS),
     1                        ((BETA(J,K),J=1,IPOS),K=1,IPOS),
     2                         (A(J),J=1,IPOS)
            READ (LUNIN) IPOS,(ALPHAQ(J),J=1,IPOS),
     1                        ((BETAQ(J,K),J=1,IPOS),K=1,IPOS),
     2                         (AQ(J),J=1,IPOS)
C
C           determining frequently used constants
C
            MP1   = M + 1
            MM2   = M * 2
            MM2P1 = MM2 + 1
            MM2M1 = MM2 - 1
            MM2M2 = MM2 - 2
         ENDIF
C
C*****************************************************************
C*            P E R F O R M I N G  O N E  S T E P                *
C*****************************************************************
C
C        the following REPEAT-loop is executed until the required
C        accuracy is achieved, i.e., until the step was executed
C        successfully. Exception: step width H falls below EPSERR,
C        which indicates that the procedure does not converge
C
 3000    CONTINUE
            IF (ORDCH) THEN
C
C              the order is increased by 1. ALPHA, BETA and A
C              are overwritten with the coefficients ALPHAQ,
C              BETAQ and AQ
C
               ORDCH = .FALSE.
               M     = M + 1
               MP1   = M + 1
               DO 160 I = 1,M
                  ALPHA(I) = ALPHAQ(I)
                  DO 150 J = 1,M
                     BETA(J,I) = BETAQ(J,I)
  150             CONTINUE
                  A(I) = AQ(I)
  160          CONTINUE
C
C              the new coefficients for ALPHAQ, BETAQ and AQ
C              are read in
C
               READ (LUNIN) IPOS,(ALPHAQ(J),J=1,IPOS),
     1                       ((BETAQ(J,K),J=1,IPOS),K=1,IPOS),
     2                        (AQ(J),J=1,IPOS)
C
C              determining frequently used constants
C
               MM2   = M * 2
               MM2P1 = MM2 + 1
               MM2M1 = MM2 - 1
               MM2M2 = MM2 - 2
               IF (.NOT. ZERO) THEN
C
C                 DELTAK and HEPS have to be determined for the
C                 new order. The old values may be used in part.
C
                  DELTAK(MM2M1) = IRKRT (N,DFDY,DFDX,WGTSUM,G,F1)
                  DELTAK(MM2)   = IRKRT (N,DFDY,DFDX,WGTSUM,G,F1)
                  IF (DELTAK(MM2M2) .NE. 0.0D0) THEN
                     HEPS(M) = (EPS * FAC(M) / DELTAK(MM2M2))
     1                          ** (1.0D0 / DBLE(MM2))
                     H       = SG * SF * HEPS(M)
                  ELSE
                     ZERO = .TRUE.
                  ENDIF
               ENDIF
               IF (ZERO) THEN
C
C                 When determining the step size HEPS(M) we have
C                 encountered division by zero. In this case we set
C                 DELTAK(L) = EPS * 10**(5-L), for L = 1, ..., 2*M-1
C                 and chose 0.1 as the initial step size.
C
                  TPOW = 1.0D0
                  DO 170 I = 1,MM2
                     DELTAK(MM2-I) = EPS * TPOW
                     TPOW          = TPOW * 10.0D0
  170             CONTINUE
                  H = SG * 0.1D0
               ENDIF
            ENDIF
C
C*****************************************************************
C*                      i t e r a t i o n                        *
C*****************************************************************
C
C           Initializing the weights KI and KIQ, i.e., DBLEK and
C           DBLEKQ to be equal to the value H*F0 for the iteration
C
            DO 190 I = 1,N
               DO 180 J = 1,M
                  DBLEK(I,J)  = H*F0(I)
                  DBLEKQ(I,J) = H*F0(I)
  180          CONTINUE
               DBLEKQ(I,MP1) = H*F0(I)
  190       CONTINUE
C
C           store the old approximation,
C           in order to be able to find an error
C           estimate after the first iteration
C
            DO 200 I = 1,N
               YOLD(I) = Y0(I) + H*F0(I)
  200       CONTINUE
C
C           In the following REPEAT-loop we perform a complete
C           iteration until the required accuracy is achieved
C           or until the number of iterations exceeds 2*M+1.
C           If the latter is the case, the order is increased
C           and the step is repeated. If a convergence of this
C           iteration is not occuring, or if the order cannot
C           be increased further, the whole step is repeated
C           with a lower step size
C
            L = 0
            FACLP1 = 1.0D0
 4000       CONTINUE
C
C              check number of function calls
C
               IF (IFU .GT. MXCALL) THEN
C
C                 if necessary stop and set error number
C
                  STOP = .TRUE.
                  IERR = 3
               ENDIF
               L      = L + 1
               FACLP1 = FACLP1 * (L+1)
               ICAUSE = 0
C
C*****************************************************************
C*           L-th full iteration to find the weights DBLEK     *
C*****************************************************************
C
               DO 240 K = 1,M
C
C                 determine intermediate points in the
C                 interval [X0, X0+H]
C
                  X1 = X0 + H*ALPHA(K)
                  DO 220 I = 1,N
                     Y(I) = Y0(I)
                     DO 210 J = 1,M
                        Y(I) = Y(I) + DBLEK(I,J)*BETA(J,K)
  210                CONTINUE
  220             CONTINUE
C
C                 insert the intermediate points into the
C                 right-hand side of the differential equation
C
                  CALL DES (N,X1,Y,F1)
C
C                 perform the full step
C
                  DO 230 I = 1,N
                     DB(I,K) = H*F1(I)
  230             CONTINUE
  240          CONTINUE
C
C              store the new weights DBLEK after the
C              full step
C
               DO 260 I = 1,N
                  DO 250 J = 1,M
                     DBLEK(I,J) = DB(I,J)
  250             CONTINUE
  260          CONTINUE
C
C*****************************************************************
C*           L-th iteration for the weights DBLEKQ using a       *
C*           complete step                                       *
C*****************************************************************
C
               DO 300 K = 1,MP1
C
C                 determine the intermediate points in the
C                 interval [X0, X0+H]
C
                  X1 = X0 + H*ALPHAQ(K)
                  DO 280 I = 1,N
                     YQ(I) = Y0(I)
                     DO 270 J = 1,MP1
                        YQ(I) = YQ(I) + DBLEKQ(I,J)*BETAQ(J,K)
  270                CONTINUE
  280             CONTINUE
C
C                 insert the intermediate points into the
C                 right-hand side of the differential equation
C
                  CALL DES (N,X1,YQ,F1)
C
C                 perform one full step
C
                  DO 290 I = 1,N
                     DB(I,K) = H*F1(I)
  290             CONTINUE
  300          CONTINUE
C
C              store the new weights DBLEK after the
C              complete step
C
               DO 320 I = 1,N
                  DO 310 J = 1,MP1
                     DBLEKQ(I,J) = DB(I,J)
  310             CONTINUE
  320          CONTINUE
C
               IFU = IFU + MM2P1
C
C*****************************************************************
C*           approximations  from  the  L-th  iteration          *
C*****************************************************************
C
C              determine new approximations Y with
C              order M, and YQ with order M+1, in the L-th
C              iteration using RUNGE-KUTTA
C
               DO 340 I = 1,N
                  Y(I)  = Y0(I)
                  YQ(I) = Y0(I)
                  DO 330 J = 1,M
                     Y(I)  = Y(I) + A(J)  * DBLEK(I,J)
                     YQ(I) = YQ(I)+ AQ(J) * DBLEKQ(I,J)
  330             CONTINUE
                  YQ(I) = YQ(I) + AQ(MP1) * DBLEKQ(I,MP1)
  340          CONTINUE
C
C              determine the absolute and relative average
C              difference DELTAG and DGREL between two
C              successive iterations, as well as the absolute
C              or relative average difference E or EREL
C              between the two approximative solutions
C              obtained using methods of different order
C
               DELTAG = 0.0D0
               E      = 0.0D0
               YQNORM = 0.0D0
               DO 350 I = 1,N
                  YQNORM  = YQNORM + YQ(I)*YQ(I)
                  YDIFF   = YQ(I) - YOLD(I)
                  DELTAG  = DELTAG + G(I) * YDIFF*YDIFF
                  YDIFF   = Y(I) - YQ(I)
                  E       = E + G(I) * YDIFF*YDIFF
                  YOLD(I) = YQ(I)
  350          CONTINUE
               DELTAG = DSQRT(DELTAG/WGTSUM)
               E      = DSQRT(E/WGTSUM)
               IF (YQNORM .GT. 0.0D0) THEN
                   YQNORM = DSQRT(YQNORM)
                   DGREL  = DELTAG/YQNORM
                   EREL   = E/YQNORM
               ELSE
                   DGREL = DELTAG
                   EREL  = E
               ENDIF
               ERROR  = DMAX1(EREL,DGREL)
               EPSLOC = DMAX1(EPSLOC,ERROR)
C
C*****************************************************************
C*             test  for  stopping the  iteration                *
C*****************************************************************
C
               LH = L
               IF (EREL .GE. EPS) THEN
C
C                 the difference between two approximations for
C                 different orders differ by a term of the same
C                 order 2*M+1 when using full steps as that of the
C                 approximate solution.
C                 If this is not the case the step size is decreased
C                 according to the theory.
C
                  ICAUSE = 1
                  IF (IC .NE. 0) SF = SF * 0.9D0
                  IC   = IC + 1
                  HBEG = H
                  H    = H * SF * (0.5D0*EPS/E)**(1.0D0/DBLE(MM2P1))
                  L    = 0
               ELSE
C
C                 compute the estimated average iteration error
C
                  IF (ZERO) THEN
                     DK = DELTAK(L-1)
                  ELSE
                     DK = 5.0D0 * H**(L+1) * DELTAK(L-1) / FACLP1
                  ENDIF
C
C                 When computing H only the local error is considered,
C                 hence we need to test for convergence of the
C                 iteration procedure
C
                  IF (DELTAG .GE. DK) THEN
                     ICAUSE = 2
                     HBEG   = H
                     H      = 0.6D0 * H
                     SF     = 0.8D0 * SF
                     L      = 0
                  ELSE
                     IF (DGREL .GE. EPS .AND. M .GE. MMAXM1 .AND.
     1                       L .GE. MM2M1) THEN
C
C                       if the maximum number of iterations is
C                       reached but without the desired accuracy,
C                       then the order should be increased. If this
C                       is not possible, since e.g. no more nodes
C                       are available, we attempt to come to a successful
C                       end by decreasing the step size
C
                        ICAUSE = 3
                        HBEG   = H
                        H      = 0.8D0 * H
                        L      = 0
                     ENDIF
                  ENDIF
               ENDIF
               IF (L .EQ. 0) THEN
C
C                 L=0 indicates that the step  cannot be
C                 completed successfully
C
                  IF (LUNPR .GT. 0) THEN
                     WRITE (LUNPR,6200) ISC,M,HBEG,X,YQ(1),ERROR,
     1                                  IFU,LH,ICAUSE
                     WRITE (LUNPR,6300) (YQ(I),I=2,N)
                  ENDIF
                  IF (H .LT. EPSERR) THEN
C
C                    the procedure does not converge
C
                     STOP = .TRUE.
                     IERR = 1
                  ELSE
C
C                    the calculations of the last step are
C                    cancelled and the step is repeated
C
                     FACLP1 = 1.0D0
                     X      = X0 + H
                     IF (SG*(XEND-X) .LT. 0.0D0) THEN
                        H = XEND - X0
                        X = XEND
                     ENDIF
                     DO 370 I = 1,N
                        DO 360 J = 1,M
                           DBLEK(I,J)  = H * F0(I)
                           DBLEKQ(I,J) = H * F0(I)
  360                   CONTINUE
                        DBLEKQ(I,MP1) = H * F0(I)
  370                CONTINUE
                     DO 380 I = 1,N
                        YOLD(I)  = Y0(I) + H*F0(I)
  380                CONTINUE
                  ENDIF
               ENDIF
C
C*****************************************************************
C*                l o o p    i n q u i r i e s                   *
C*****************************************************************
C
            IF (((L .LT. MM2M1 .AND. DGREL .GE. EPS) .OR. L .EQ. 0)
     1                         .AND. .NOT. STOP) GOTO 4000
            IF (DGREL .GE. EPS .AND. MP1 .LT. MMAX
     1                         .AND. .NOT. STOP) THEN
C
C              the desired accuracy has not been reached after
C              the theoretically maximal number of iterations. The
C              order is increased and the step is repeated
C
               ORDCH  = .TRUE.
               ICAUSE = 4
               IF (LUNPR .GT. 0) THEN
                  WRITE (LUNPR,6200) ISC,M,H,X,YQ(1),ERROR,IFU,
     1                               L,ICAUSE
                  WRITE (LUNPR,6300) (YQ(I),I=2,N)
               ENDIF
            ELSE
C
C              prevent from overstep the maximal order
C
               ORDCH = .FALSE.
            ENDIF
         IF (ORDCH) GOTO 3000
         IF (.NOT. STOP) THEN
C
C           the step was successfully completed
C
            IF (LUNPR .GT. 0) THEN
               WRITE (LUNPR,6200) ISC,M,H,X,YQ(1),ERROR,IFU,L,
     1                            ICAUSE
               WRITE (LUNPR,6300) (YQ(I),I=2,N)
            ENDIF
            IF (LUNOUT .GT. 0) THEN
                  I = 1
                  WRITE (LUNOUT,7200) ISC,X,I,YQ(1),ERROR
                  WRITE (LUNOUT,7300) (I,YQ(I),I=2,N)
            ENDIF
         ENDIF
C
         IF (X .NE. XEND .AND. .NOT. STOP) THEN
C
C           the end of the integration interval has not
C           been reached
C
            IF (IC .GT. 1) SF = SF / 0.97D0
            X0 = X
            DO 390 I = 1,N
               Y0(I) = YQ(I)
  390       CONTINUE
         ELSE
C
C           program stop since the initial value problem was
C           solved (IERR=0), or because the procedure does not
C           converge (IERR=1)
C
            STOP  = .TRUE.
            EPS   = EPSLOC
         ENDIF
      IF (.NOT. STOP) GOTO 1000
C
C*****************************************************************
C*           f o r m a t  s t a t e m e n t s                    *
C*****************************************************************
C
 5100 FORMAT ('1','INITIAL CONDITION:',//,
     1        1X,10X,'X0',10X,2X,'COMP.',2X,10X,'Y0',/,
     2        1X,E22.15,      1X,I4,1X  ,2X,E22.15)
 5200 FORMAT (1X,            23X,I4,1X  ,2X,E22.15)
 5300 FORMAT ('0','RIGHT ENDPOINT OF THE INTERVAL OF INTEGRATION: ',
     1        E22.15,/,30X,'DESIRED ACCURACY: ',E22.15,/
     2         ,25X,'LARGEST ORDER ALLOWED:',I4)
 5350 FORMAT ('0','REASON FOR DECREASE IN STEP SIZE' ,/,5X
     2        ,'0  NO DECREASE',/,5X
     3        ,'1  EREL >= EPS',/,5X,'2  DELTAG >= DK',/,5X
     4        ,'3  DGREL >= EPS')
 5400 FORMAT ('0','WEIGHTS G:       COMP.  ',10X,'G',/,
     1         (17X,I4,3X,E22.15))
 6000 FORMAT ('0','NUMBER OF THE OUTPUT FILE (0=NONE):',I4)
 6100 FORMAT ('0',' STEP  ORDER STEP SIZE  UPPER BOUND APPROXIMATION',
     1        ' ERROR      FUNCTION ITERATION CAUSE',/,18X,'H',24X,
     2        'Y',7X,'ESTIMATE',3X,'CALLS',4X,'STEPS')
 6200 FORMAT (1X,I6,2X,I3,2X,E10.3,1X,E10.3,2X,E13.6,1X,E10.3,1X,
     1        I8,1X,I9,3X,I1)
 6300 FORMAT (36X,E13.6)
 7000 FORMAT ('0','NUMBER OF THE LOG FILE (0=NONE):',I4)
 7100 FORMAT ('0',' STEP  BOUND FOR UPPER LIMIT  COMP.',2X,
     1        'APPROXIMATE SOLUTION     ERROR ESTIMATE',/,7X,
     2        ' OF INTEGRATIONINTERVAL ')
 7200 FORMAT (1X,I6,1X,E22.15,1X,I4,2X,D24.17,1X,E18.11)
 7300 FORMAT (31X,I4,2X,D24.17)
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION IRKRT (N,DFDY,DFDX,WGTSUM,G,DFDXH)
C[IX{IRKRT}*)
C
C*****************************************************************
C                                                                *
C The DOUBLE PRECISION FUNCTION IRKRT determines square roots,   *
C that occur in SUBROUTINE IMRUKU when determining the step size *
C H and the estimated average iteration error DK or DELTAK(L),   *
C L=0, 1, ...., 2*M+1.                                           *
C[BE*)
C                                                                *
C IRKRT = SQRT(G*(DF/DX*(DF/DY)**K))**2/WGTSUM), for K=1,..,2*M-1*
C                                                                *
C Here we determine (DF/DY)**K during each call of IRKRT, by     *
C storing (DF/DY)**(K-1)*DF/DX in vector previously used for     *
C DF/DX and DFDX.                                                *
C                                                                *
C                                                                *
C INPUT PARAMETERS: (compare description of variables in IMRUKU) *
C =================                                              *
C N       : dimension                                            *
C DFDY    : 2-dimensional array DFDY(1:N,1:N); the derivative    *
C           DF/DY of the right-hand side of the differential     *
C           equation with respect to Y from IMRUKU at (X0,Y0).   *
C DFDX    : vector DFDX(1:N); the derivative DF/DX of the right- *
C           hand side of the differential equation with respect  *
C           to X from IMRUKU or (DF/DY)**(K-1) * DF/DX at (X0,Y0)*
C G       : vector G(1:N) of weights                             *
C WGTSUM  : G(1) + G(2) + ... + G(N)                             *
C DFDXH   : vector DFDXH(1:N); auxiliary storage                 *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C IRKRT   : computed square root (compare program description)   *
C DFDX    : vector DFDX(1:N);  (DF/DY)**K * DF/DX                *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C I,J     : control variables                                    *
C SUM     : scalar products DFDY(I,..)*DFDX                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Thomas Eul                                         *
C  date     : 08.12.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     parameter
C
      INTEGER N
      DIMENSION DFDX(N), DFDY(N,N), G(N), DFDXH(N)
      DOUBLE PRECISION  WGTSUM
C
C     local variables
C
      INTEGER I, J
      DOUBLE PRECISION  SUM
C
      DO 20 I = 1,N
         SUM = 0.0D0
         DO 10 J = 1,N
            SUM = SUM + DFDY(I,J) * DFDX(J)
   10    CONTINUE
         DFDXH(I) = SUM
   20 CONTINUE
      SUM = 0.0D0
      DO 30 J = 1,N
         DFDX(J) = DFDXH(J)
         SUM = SUM + DFDX(J)**2 * G(J)
   30 CONTINUE
      IRKRT = DSQRT(SUM/WGTSUM)
C
      RETURN
      END

![  {Implicit Runge--Kutta Methods of Gaussian Type}                    
![  {Implicit Runge--Kutta Methods of Gaussian Type}*)                  
      SUBROUTINE IRKDRV (DES, N, MMAX, IFLAG, LUNIN, LUNOUT, LUNPR,     &
      IERR, EPSM, EPS, G, X0, XEND, Y0, YQ, WORK1, WORK21, WORK22,      &
      WORK23)                                                           
!                                                                       
!*****************************************************************      
!                                                                *      
! IRKDRV solves an initial value problem (IVP) with N            *      
! differential equations of 1st order.                           *      
! The IVP is solved by the implicit RUNGE-KUTTA method (IRKV).   *      
! In it both step size control and a control of the order of the *      
! IRKV is used. The nodes for the IRKV are stored in a file,     *      
! numbered LUNIN. This file can be produced by the SUBROUTINE    *      
! IRKCOE.                                                        *      
! SUBROUTINE IRKDRV is only a driver routine for providing       *      
! storage space, for testing of the input parameters and for     *      
! calling two subroutines according to the setting of the        *      
! parameter IFLAG:                                               *      
! IFLAG=0       : We only solve the SVP. We assume that the node *      
!                 file already exists.                           *      
! IFLAG=1       : Before solving the IVP we must produce the node*      
!                 file numbered LUNIN, by SUBROUTINE IRKCOE.     *      
!                                                                *      
! The solution of the IVP is found in the SUBROUTINE IMRUKU.     *      
! A description of the algorithm used can be found there.        *      
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
! IFLAG   : =0; solve the IVP only.                              *      
!           =1; before solving the IVP, we must produce the node *      
!               file. The generated nodes are the GAUSS-LEGENDRE *      
!               nodes.                                           *      
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
!           If EPSM <= 0, the machine constant is taken from     *      
!           IRKDVR.                                              *      
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
! WORK1   : auxiliary vector WORK1(1:8*MMAX+5*N-2)               *      
! WORK21  : 2-dimensional auxiliary array WORK21(1:N,1:N)        *      
! WORK22  : 2-dimensional auxiliary array                        *      
!           WORK22(1:N+MMAX-1,1:MMAX-1)                          *      
! WORK23  : 2-dimensional auxiliary array                        *      
!           WORK23(1:2*N+MMAX,1:MMAX)                            *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! IERR    : error parameter                                      *      
!           =0; run was successful                               *      
!           =1; no convergence. Possible remedy: increase the    *      
!               maximal order MMAX.                              *      
!           =2; wrong input parameter(s)                         *      
! EPS     : estimate for the largest local relative error        *      
! YQ      : vector YQ(1:N); approximation for the solution of    *      
!           the IVP.                                             *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I       : loop variable                                        *      
! FMACHP  : auxiliary variable for determining the machine       *      
!           constant                                             *      
! ZERO    : logic variable                                       *      
!           =.TRUE. ,if G(I)=0 for all I                         *      
!           =.FALSE. otherwise                                   *      
!                                                                *      
!                                                                *      
! All other local variables are INTEGER variables and serve for  *      
! subdividing the auxiliary arrays for the subroutines:          *      
!                                                                *      
! IFAC, IDELTK, IHEPS, IALPH, IALPHQ, IA, IAQ, IF0, IF1, IDFDX,  *      
! IY, IYOLD, IDBLK, IBETA, IBETAQ, IDBLKQ, IDB, IC.              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: IRKCOE, IMRUKU, IRKRT, GALE0, MACHPD    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 09.17.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!     parameters                                                        
!                                                                       
      EXTERNAL DES 
      INTEGER N, MMAX, IFLAG, LUNIN, LUNOUT, LUNPR, IERR 
      DIMENSION G (N), Y0 (N), YQ (N), WORK1 (8 * MMAX + 5 * N - 2),    &
      WORK21 (N, N), WORK22 (N + MMAX - 1, MMAX - 1), WORK23 (2 * N +   &
      MMAX, MMAX)                                                       
      DOUBLEPRECISION EPSM, EPS, X0, XEND 
!                                                                       
!     local variables                                                   
!                                                                       
      LOGICAL ZERO 
      INTEGER I, IFAC, IDELTK, IHEPS, IALPH, IALPHQ, IA, IAQ, IF0, IF1, &
      IDFDX, IY, IYOLD, IDBLK, IBETA, IBETAQ, IDBLKQ, IDB, IC, IHELP,   &
      IBETA2, IDBKQ2, IDB2                                              
      DOUBLEPRECISION FMACHP 
!                                                                       
!*****************************************************************      
!*               testing the input parameters                    *      
!*****************************************************************      
!                                                                       
      I = 1 
      ZERO = .TRUE. 
   10 CONTINUE 
      IF (G (I) .NE.0.0D0) ZERO = .FALSE. 
      I = I + 1 
      IF (ZERO.AND.I.LE.N) GOTO 10 
!                                                                       
      IF (N.LE.0.OR.MMAX.LE.4.OR.LUNIN.LE.0.OR.LUNOUT.LT.0.OR.LUNPR.LT.0&
     &.OR.EPS.LE.0.0D0.OR.ZERO.OR.IFLAG.LT.0.OR.IFLAG.GT.1) THEN        
         IERR = 2 
      ELSE 
!                                                                       
!*****************************************************************      
!*  determine the addresses for subdividing the auxiliary arrays *      
!*****************************************************************      
!                                                                       
         IFAC = 1 
         IDELTK = IFAC + MMAX 
         IHEPS = IDELTK + 2 * MMAX 
         IALPH = IHEPS + MMAX 
         IALPHQ = IALPH + MMAX - 1 
         IA = IALPHQ + MMAX 
         IAQ = IA + MMAX - 1 
         IF0 = IAQ + MMAX 
         IF1 = IF0 + N 
         IDFDX = IF1 + N 
         IY = IDFDX + N 
         IYOLD = IY + N 
!                                                                       
         IDBLK = 1 
         IBETA = IDBLK + N * (MMAX - 1) 
         IHELP = N + MMAX - 1 
         IF (IHELP.LT.IBETA) THEN 
            IBETA2 = INT (IBETA / IHELP) 
            IBETA = MOD (IBETA, IHELP) 
            IF (IBETA.EQ.0) THEN 
               IBETA = IHELP 
            ELSE 
               IBETA2 = IBETA2 + 1 
            ENDIF 
         ELSE 
            IBETA2 = 1 
         ENDIF 
         IBETAQ = 1 
         IDBLKQ = IBETAQ + MMAX * MMAX 
         IHELP = 2 * N + MMAX 
         IF (IHELP.LT.IDBLKQ) THEN 
            IDBKQ2 = INT (IDBLKQ / IHELP) 
            IDBLKQ = MOD (IDBLKQ, IHELP) 
            IF (IDBLKQ.EQ.0) THEN 
               IDBLKQ = IHELP 
            ELSE 
               IDBKQ2 = IDBKQ2 + 1 
            ENDIF 
         ELSE 
            IDBKQ2 = 1 
         ENDIF 
         IDB = IBETAQ + MMAX * (MMAX + N) 
         IF (IHELP.LT.IDB) THEN 
            IDB2 = INT (IDB / IHELP) 
            IDB = MOD (IDB, IHELP) 
            IF (IDB.EQ.0) THEN 
               IDB = IHELP 
            ELSE 
               IDB2 = IDB2 + 1 
            ENDIF 
         ELSE 
            IDB2 = 1 
         ENDIF 
!                                                                       
!                                                                       
!*****************************************************************      
!*       create the node file, if IFLAG > 0                      *      
!*****************************************************************      
!                                                                       
         IF (IFLAG.GT.0) THEN 
            IC = 1 
!                                                                       
            CALL IRKCOE (MMAX, LUNIN, WORK1 (IC), WORK1 (IAQ), WORK1 (  &
            IALPHQ), WORK23 (IBETAQ, 1) )                               
         ENDIF 
         IF (EPSM.LE.0.0D0) THEN 
!                                                                       
!*****************************************************************      
!*             determine the machine constant                    *      
!*****************************************************************      
!                                                                       
!           EPSM ist the smallest positive machine number with          
!           (1.0+EPSM) .GT. 1.0 . EPSM  is determined                   
!           approximately as a power of 1./2. .                         
!                                                                       
            FMACHP = 1.0D0 
   20       CONTINUE 
            FMACHP = 0.5D0 * FMACHP 
            IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 20 
            EPSM = 2.0D0 * FMACHP 
         ENDIF 
!                                                                       
!*****************************************************************      
!*              solve the initialt value problem                 *      
!*****************************************************************      
!                                                                       
         CALL IMRUKU (DES, N, MMAX, LUNIN, LUNOUT, LUNPR, IERR, EPSM,   &
         EPS, G, X0, XEND, Y0, YQ, WORK1 (IFAC), WORK1 (IDELTK),        &
         WORK1 (IHEPS), WORK1 (IALPH), WORK1 (IALPHQ), WORK1 (IA),      &
         WORK1 (IAQ), WORK1 (IF0), WORK1 (IF1), WORK1 (IDFDX), WORK1 (  &
         IY), WORK1 (IYOLD), WORK21, WORK22 (IDBLK, 1), WORK22 (IBETA,  &
         IBETA2), WORK23 (IBETAQ, 1), WORK23 (IDBLKQ, IDBKQ2), WORK23 ( &
         IDB, IDB2) )                                                   
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE IRKDRV                         

C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 15.12}
C[  {Adaptive Quadrature Methods}
C[  {Adaptive Quadrature Methods}*)
C[LE*)
      SUBROUTINE GAX (INTVAL,EPS,N,FCT,NMAX,LQM,XNODES,
C[IX{GAX}*)
     +                QVAL,EXE,TPOINT,STATUS,IERR,AK,IDGR)
C
C*****************************************************************
C                                                                *
C            Adaptive quadrature in one dimension.               *
C[BE*)
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  INTVAL: 2-dimensional array INTVAL(1:100,1:2); user supplied  *
C          sub-intervals of integration:                         *
C                        INTVAL(1,1) = A1, INTVAL(1,2) = B1      *
C                           ...              ...                 *
C                        INTVAL(N,1) = AN, INTVAL(N,2) = BN.     *
C          Here AI must be less than BI for each I.              *
C                                                                *
C  EPS   : relative accuracy required for the solution           *
C                                                                *
C  N     : number of initially supplied sub-intervals            *
C                                                                *
C  FCT   : name of the function FCT(X) that is to be integrated, *
C          provided by the user in the form                      *
C               DOUBLE PRECISION FUNCTION FCT (X).               *
C          In the calling program it must be defined as EXTERNAL.*
C                                                                *
C  NMAX  : maximal number of sub-intervals allowed to be produced*
C          during the calculations                               *
C                                                                *
C  LQM   : parameter that determines which quadrature formula    *
C          to use on the sub-intervals                           *
C             LQM = 1: trapezoidal rule                          *
C                 = 2: GAUSS quadrature formula                  *
C                 = 3: CLENSHAW-CURTIS formula                   *
C                 = 4: ROMBERG method                            *
C                 = 5: summed NEWTON-COTES formulas              *
C                                                                *
C  IDGR  : this label determines the degree of the quadrature    *
C          formula used                                          *
C          LQM =1:IDGR not used                                  *
C          LQM =2:IDGR = 2: GAUSS quadrature formula of degree 2 *
C                 IDGR = 3: GAUSS quadrature formula of degree 3 *
C                 ...                                            *
C                 IDGR =20: GAUSS quadrature formula of degree 20*
C                 ( 2 <= IDGR <= 20 )                            *
C          LQM =3:IDGR = 2: CLENSHAW-CURTIS formula with         *
C                           2 + 1 weights                        *
C                 IDGR = 4: CLENSHAW-CURTIS formula with         *
C                           4 + 1 weights                        *
C                 ...                                            *
C                 IDGR must be greater than 1 and even.          *
C          LQM =4:IDGR not used                                  *
C          LQM =5:IDGR = 2: summed NEWTON-COTES formula for 2    *
C                           nodes, i.e., SIMPSON's rule          *
C                 IDGR = 3: summed NEWTON-COTES formula for 3    *
C                           nodes, i.e., the 3/8 formula         *
C                 ...                                            *
C                 IDGR = 7: summed NEWTON-COTES formula for 7    *
C                           nodes, i.e., the 7/17280 formula     *
C                 (We must have 2 <= IDGR <= 7 here)             *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C                                                                *
C  XNODES: vector XNODES(1:NMAX * 4) containing the sub-intervals*
C          produced during the calculations in the following     *
C          form:                                                 *
C          XNODES(1) : error estimate for the quadrature from A  *
C                      to B                                      *
C          XNODES(2) : computed value of the integral from A to B*
C          XNODES(3) : left hand endpoint A of the interval      *
C          XNODES(4) : right hand endpoint B of the interval     *
C                                                                *
C  QVAL  : approximate value for the integral when the procedure *
C          is stopped.                                           *
C                                                                *
C  EXE   : accuracy of QVAL achieved when stopping               *
C          (absolute error)                                      *
C                                                                *
C  NMAX  : actual number of sub-intervals used                   *
C                                                                *
C  IERR  : error parameter:                                      *
C          = 0: everything is o.k.,  EXE < EPS                   *
C          = 1: exceeding  MAXLEV = 100 while                    *
C               the required accuracy could not be achieved      *
C               => choose a larger value for IDGR or a smaller   *
C                  accuracy bound EPS                            *
C          = 2: N > NMAX, the required accuracy was not achieved *
C          = 3: left endpoint A > right endpoint B:              *
C               endpoints entered incorrectly have been exchanged*
C          = 4: conditions IERR = 3 and = 1 both hold            *
C          = 5: conditions IERR = 3 and = 2 both hold            *
C          = 6: degree IDGR of the GAUSS quadrature formula      *
C               satisfies 2 > IDGR > 20                          *
C          = 7: weights of the CLENSHAW-CURTIS formula           *
C               IDGR are odd or less than 2                      *
C          = 8: number of nodes IDGR in the NEWTON-COTES formula *
C               2 > IDGR > 7                                     *
C          = 9: N < 0  or  EPS < 1.D-12  or  NMAX < 1            *
C          =10: 1 > LQM > 5                                      *
C          ....                                                  *
C                                                                *
C          >=100: error when determining an integral value       *
C                                                                *
C  AUXILIARY VECTORS:                                            *
C  ==================                                            *
C  ITNODE: vector ITNODE(1:2); contains information on the       *
C          current sub-interval:                                 *
C          ITNODE(1) indicates the state of the sub-interval:    *
C                    = 1: state = small                          *
C                    = 0: state = big                            *
C          ITNODE(2) contains the current level of the sub-      *
C                    interval                                    *
C  EL    : vector EL(1:102); contains the current row of the     *
C          ROMBERG scheme for RICHARDSON extrapolation           *
C  ELEPS : vector ELEPS(1:102); contains the current row of the  *
C          EPSILON-algorithm                                     *
C  TNODE : vector TNODE(1:4); contains the sub-interval currently*
C          in use                                                *
C  XZERO : vector XZERO(1:21); contains the zeros of the LEGENDRE*
C          polynomials in ascending order                        *
C  ZWGH  : vector ZWGH(1:21); contains the weights corresponding *
C          to the zeros of the GAUSS-quadrature formulas         *
C                                                                *
C  AK    : 2-dimensional array AK(0:IDGR, 2); contains the       *
C          weights and nodes of the CLENSHAW-CURTIS quadrature   *
C          formula for the reference interval [-1,1] that were   *
C          produced during the calculations                      *
C  TPOINT: vector TPOINT(1:NMAX); used as a marker:              *
C          this vector contains XNODES in quasi order; TPOINT(1) *
C          points to the sub-interval with the largest quadrature*
C          error, etc.                                           *
C  STATUS: vector STATUS(1:NMAX * 2); contains the state of each *
C          sub-interval (see vector ITNODE(1:4) )                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required :                                        *
C       GXENT   : tests the input parameters                     *
C       GXDIV   : splits the intervals                           *
C       GXQUAD  : determines the integral approximately:         *
C          GXGAUS  : integration using GAUSS quadrature          *
C          CLENSH  : integration using CLENSHAW-CURTIS formulas  *
C          QUAROM  : integration using ROMBERG method            *
C          QUANEC  : integration using NEWTON-COTES formulas     *
C             FCT  : function to be integrated                   *
C       GXINS   : inserts new sub-intervals                      *
C       GXDEL   : erases a sub-interval                          *
C       GXACC   : accesses a sub-interval                        *
C       GXRIEP  : RICHARDSON-extrapolation, EPSILON-algorithm    *
C       GALE0   : nodes and weights for GAUSS-quadrature         *
C          GXPOLY  : HORNER scheme                               *
C          GXPEGA  : PEGASUS method                              *
C       WGKNOT  : weights, nodes for CLENSHAW-CURTIS formulas    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  sources : 1. Eul, T. and  Rheinbach, H.J.:                    *
C               Lecture: The combination of adaptive and extra-  *
C               polation methods for numerical integration,      *
C               RWTH Aachen, October 1984.                       *
C            2. D. Kahaner and J. Stoer, see [KAHA83].           *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      EXTERNAL         FCT
      INTEGER          IDGR, GRDFLG, IERR2
      DOUBLE PRECISION TNODE(4), N1(4), N2(4), XNODES(1:NMAX, 1:4)
      DOUBLE PRECISION EL(0:102), ELEPS(0:102), INTVAL(1:100,1:2)
      DOUBLE PRECISION XZERO(21), ZWGH(21),AK(0:IDGR, 2), HELP, EPS
      INTEGER          TPOINT(NMAX), STATUS(1:NMAX, 1:2), ITNODE(2)
      DOUBLE PRECISION QVAL, EXE, QEXE, EPSL, EPSEXE, BEE, TEE
C
C     testing the input parameters,
C     if necessary return
C
      IERR = 0
      CALL GXENT (LQM, N, NMAX, EPS, AK, IDGR, ZWGH, XZERO,
     +            GRDFLG, IERR)
      IF (IERR .GT. 0) RETURN
C
C     initialize the pointer vector TPOINT
C
      DO 10 I = 1, NMAX
         TPOINT(I) = I
 10   CONTINUE
C
C     initialize the vector XNODES for future calculations
C
      QVAL     = 0.0D0
      TEE      = 0.0D0
      TNODE(1) = 1000.0D0
      DO 30 I = 1, N
C
C        the integral is computed for the I-th
C        starting interval
C
         TNODE(3) = INTVAL(I,1)
         TNODE(4) = INTVAL(I,2)
C
C        check that left endpoint < right endpoint
C
         IF (DABS(TNODE(3)) .GT. DABS(TNODE(4))) THEN
C
C           if right endpoint < left endpoint:
C           ==> exchange the two
C
            IERR     = 3
            HELP     = TNODE(3)
            TNODE(3) = TNODE(4)
            TNODE(4) = HELP
         ENDIF
C
         CALL GXQUAD (TNODE, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR,
     +                IERR2)
C
C        store the I-th sub-interval in the vector XNODES
C
         DO 40 J = 1, 4
            XNODES (I,J) = TNODE(J)
 40      CONTINUE
C
         STATUS(I,1) = 1
         STATUS(I,2) = 0
         QVAL        = QVAL + TNODE(2)
         TEE         = TEE + 1000.0D0
 30   CONTINUE
C
C     initialize remaining variables
C
      LEVEL    = 0
      MAXLEV   = 100
      BEE      = 0.0D0
      EXE      = TEE
      EPSEXE   = TEE
      EL(1)    = QVAL
      ELEPS(1) = QVAL
      EPSL     = 0.0D0
      QEXE     = 1000.0D0
      M        = 1
C
C     integrate
C
      DO 111 ITER = 1, NMAX
C
C        test break-off criterion of the method used
C
         IF ( TEE .LE. EPSL .OR. EXE .LE. EPSL .OR.
     +        LEVEL .GT. MAXLEV ) THEN
C
C           If TEE <= EPSL  or  EXE <= EPSL  or  LEVEL > MAXLEV
C           ==> correct computation
C
            IF (IERR2 .GT. 0)      IERR = IERR + IERR2
            IF (LEVEL .GT. MAXLEV) IERR = IERR + 1
            IF (EXE .LT. TEE) THEN
               QVAL = QEXE
            ELSEIF (TEE .LT. EXE) THEN
               EXE  = TEE
            ENDIF
            NMAX = N
            R E T U R N
         ENDIF
C
C        retrieve the interval with the largest error and transfer
C        it to the vector TNODE
C
         CALL GXACC (NMAX, XNODES, STATUS, N, TPOINT, TNODE, ITNODE, 1)
         CALL GXDEL (NMAX, N, TPOINT, 1)
         CALL GXDIV (TNODE, N1, N2, FCT, XZERO, ZWGH, LQM, GRDFLG,
     +               AK, IDGR, IERR2)
C
C
         IF (ITNODE(2) .EQ. LEVEL) THEN
            LEVEL = LEVEL + 1
C
C           all used sub-intervals are assigned the state 'BIG'.
C           This requires to recalculate the value for BEE
C
            BEE = 0.0D0
            DO 100 I = 1, N
               STATUS (TPOINT(I), 1) = 0
               BEE                   = BEE + XNODES (TPOINT(I), 1)
 100        CONTINUE
            ITNODE(1) = 1
            ITNODE(2) = LEVEL
         ELSE
            IF (ITNODE(2) .EQ. (LEVEL-1)) THEN
               ITNODE(1) = 1
               BEE       = DABS(BEE - TNODE(1))
            ELSE
               ITNODE(1) = 0
               BEE       = DABS(BEE + (N1(1) + N2(1) - TNODE(1)) )
            ENDIF
            ITNODE(2) = ITNODE(2) + 1
         ENDIF
         IPOS = 1
C
C        check whether two additional sub-intervals can be inserted
C
 2       IF ((N+2) .GT. NMAX) THEN
C
C           If N+2 > NMAX
C           ==>  correct computations
C
            IERR = IERR + 2
            IF (IERR2 .GT. 0) IERR = IERR + IERR2
            IF (EXE .LT. TEE) THEN
               QVAL = QEXE
            ELSEIF (TEE .LT. EXE) THEN
               EXE  = TEE
            ENDIF
            NMAX = N
            R E T U R N
         ENDIF
         CALL GXINS (NMAX, XNODES, STATUS, N, TPOINT, N1, N2, ITNODE,
     +               IPOS)
C
C        updating of values  for TEE  and  QVAL
C
         TEE  = DABS(TEE + (N1(1)+N2(1)-TNODE(1)) )
         QVAL = QVAL + (N1(2)+N2(2)-TNODE(2))
C
         IF (BEE .GT. EPSL .AND. TEE .GT. EPSL) THEN
C
C           find the sub-interval with state 'BIG'
C           that has the largest quadrature error
C
            DO 110 I = IPOS, N
               IF (STATUS (TPOINT(I), 1) .NE. 1) THEN
                  IPOS = I
                  CALL GXACC (NMAX, XNODES, STATUS, N, TPOINT,
     +                        TNODE, ITNODE, IPOS)
                  CALL GXDEL (NMAX, N, TPOINT, IPOS)
                  CALL GXDIV (TNODE, N1, N2, FCT, XZERO, ZWGH,
     +                        LQM, GRDFLG, AK, IDGR, IERR2)
C
C                 are the new rows of the state 'BIG' or 'SMALL' ?
C
                  IF (ITNODE(2) .EQ. (LEVEL-1)) THEN
                     ITNODE(1) = 1
                     BEE       = DABS(BEE - TNODE(1))
                  ELSE
                     ITNODE(1) = 0
                     BEE       = DABS(BEE + (N1(1)+N2(1)-TNODE(1)))
                  ENDIF
                  ITNODE(2) = ITNODE(2) + 1
C
                  GOTO 2
C
               ENDIF
 110        CONTINUE
         ENDIF
         CALL GXRIEP (QVAL, EL ,ELEPS, EPSEXE, EPSL, M, EPS, EXE, QEXE)
 111  CONTINUE
      IF (IERR2 .GT. 0) IERR = IERR + IERR2
      IF (EXE .LT. TEE) THEN
         QVAL = QEXE
      ELSEIF (TEE .LT. EXE) THEN
         EXE  = TEE
      ENDIF
      NMAX = ITER
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXENT (LQM, N, NMAX, EPS, AK, IDGR, ZWGH,
C[IX{GXENT}*)
     +                  XZERO, GRDFLG, IERR)
C
C*****************************************************************
C                                                                *
C     This SUBROUTINE tests the validity of the input parameters,*
C     and, if necessary, determines the weights for a CLENSHAW-  *
C     CURTIS or GAUSS quadrature formula of degree IDGR.         *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMERTERS:                                         *
C     ==================                                         *
C     LQM       : quadrature formula to be used                  *
C     N         : number of starting sub-intervals               *
C     NMAX      : maximum number of sub-intervals that can be    *
C                 produced                                       *
C     EPS       : relative accuracy bound                        *
C     IDGR      : degree of the quadrature formula               *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     AK        : 2-dimensional array AK(0:IDGR, 2), the weights *
C                 and nodes for CLENSHAW-CURTIS                  *
C     ALPHA     : vector ALPHA(1:21), the nodes for GAUSS        *
C     ZWGH      : vector ZWGH(1:21), the weights for GAUSS       *
C     XZERO     : vector XZERO(1:21), the zeros for GAUSS        *
C     GRDFLG    : =1, if 0 is a node; =0, if zero is not a node  *
C     IERR      : error parameter                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required : GALE0, WGKNOT                       *
C                                                                *
C*****************************************************************
C                                                                *
C     author   : Norbert Vogt                                    *
C     date     : 05.18.1989                                      *
C     source   : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER            N, IDGR, IERR, GRDFLG, NMAX
      DOUBLE PRECISION   AK (0:IDGR, 2), ZWGH (21), XZERO (21)
      DOUBLE PRECISION   EPS
C
C     testing the input parameter
C
      IF (N .LT. 0 .OR. EPS .LT. 1.0D-12 .OR. NMAX .LT. 1)  THEN
C
C        If  N < 0   or  EPS < 1.D-12  or  NMAX < 1
C        ==> return to calling program
C
         IERR = 9
         R E T U R N
      ELSEIF (LQM .LT. 1 .OR. LQM .GT. 5)  THEN
C
C        If  LQM < 1  or  LQM > 5
C        ==> return to calling program
C
         IERR = 10
         R E T U R N
      ELSEIF (LQM .EQ. 2) THEN
C
C        If LQM = 2
C        ==> use GAUSS quadrature
C
         IF (IDGR .LE. 1 .OR. IDGR .GT. 20) THEN
C
C           If  IDGR <= 1  or  IDGR > 20
C           ==> return to calling program
C
            IERR = 6
            R E T U R N
         ELSE
C
C           the desired GAUSS quadrature formula
C           of degree 1 < IDGR  < 20 is generated
C
            CALL GALE0 (IDGR, .TRUE., XZERO, ZWGH)
C
C           nodes
C
            GRDFLG = 0
C
C           Is 0 a node ?
C
            IF ((MOD(IDGR,2) .EQ. 1)) GRDFLG = 1
         ENDIF
      ELSEIF (LQM .EQ. 3) THEN
C
C        If  LQM = 3
C        ==> use CLENSHAW-CURTIS formulas
C
         IF (MOD(IDGR,2) .NE. 0 .OR. IDGR .LT. 2) THEN
C
C           If  IDGR is not even or if IDGR < 2
C           ==> return to calling program
C
            IERR = 7
            R E T U R N
         ELSE
C
C           the weights and nodes for the CLENSHAW-
C           CURTIS formulas are determined
C
            CALL WGKNOT (IDGR, AK, IERR2)
         ENDIF
      ELSEIF (LQM .EQ. 5) THEN
C
C        If  LQM = 5
C        ==> use summed NEWTON-COTES formula
C
         IF (IDGR .LT . 2 .OR. IDGR .GT. 7) THEN
C
C           If  IDGR < 2 or IDGR > 7
C           ==> return to calling program
C
            IERR = 8
            R E T U R N
         ENDIF
C
C        for LQM = 1 and LQM = 4 no further
C        tests are required
C
      ENDIF
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXQUAD (NODE, FCT, XZERO, ZWGH, LQM, GRDFLG,
C[IX{GXQUAD}*)
     +                    AK, IDGR, IERR2)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE determines the approximate integral of the    *
C  FUNCTION FCT(X) over the interval (A,B) by the quadrature     *
C  formula determined by LQM as follows:                         *
C        LQM = 1:  trapezoidal rule                              *
C        LQM = 2:  IDGR =2,...,20 : I-point-formula of GAUSS     *
C        LQM = 3:  IDGR =2,4,6,...: CLENSHAW-CURTIS formula with *
C                                   IDGR + 1 weights             *
C        LQM = 4:  ROMBERG method                                *
C        LQM = 5:  IDGR =2,...,7  : NEWTON-COTES formulas with   *
C                                   IDGR nodes                   *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  NODE  : contains the index of the sub-interval over which we  *
C          integrate                                             *
C  FCT   : name of integrand                                     *
C  XZERO : nodes of the GAUSS-quadrature formula                 *
C  ZWGH  : weights for the nodes of the GAUSS-formulas           *
C  LQM   : indicates quadrature formula to be used               *
C  GRDFLG: =1, if 0 is a node; =0, if 0 is not a node            *
C  AK    : weights and nodes for the CLENSHAW-CURTIS formulas    *
C  IDGR  : degree of the quadrature formula                      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  NODE(2): computed value for the integral over the interval    *
C           [NODE(3), NODE(4)]                                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : CLENSH, QUAROM, QUANEC, GXGAUSS        *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 05.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL         FCT
      INTEGER          IDGR
      DOUBLE PRECISION NODE(4), XZERO(IDGR), ZWGH(IDGR), XL
      DOUBLE PRECISION AK(0:IDGR, 2), H, EL(10), ERREST, QV
      INTEGER          GRDFLG
C
      IF (LQM .EQ. 2 ) THEN
C
C        If  LQM = 2
C        ==> the integral is determined using the GAUSS-quadrature
C            formula of degree IDGR
C
         CALL GXGAUS (NODE(3), XZERO, ZWGH, IDGR, GRDFLG, FCT, QV)
C
      ELSEIF (LQM .EQ. 3) THEN
C
C        If  LQM = 2
C        ==> the integral is determined using the summed CLENSHAW-
C            CURTIS formula with IDGR weights
C
         CALL CLENSH (FCT, IDGR, NODE(3), 1, AK, QV, IERR)
C
      ELSEIF (LQM .EQ. 4) THEN
C
C        If LQM = 4
C        ==> the integral is determined via the ROMBERG method
C
         H     = 0.0D0
         NROWS = 10
         CALL QUAROM (NODE(3), NODE(4), 1.0D-8, NROWS, H,
     +                FCT, EL, QV, ERREST, IERR)
C
      ELSEIF (LQM .EQ. 5) THEN
C
C        If  LQM = 5
C        ==> the integral is determined by the summed
C            NEWTON-COTES formula with IDGR nodes
C
         CALL QUANEC (NODE(3), NODE(4), 1, IDGR, FCT, QV,
     +                EL(1), EL(2), IERR)
C
      ELSE
C
C        If  LQM = 1
C        ==> the integral is determined by the trapezoidal rule
C
         XL   = (NODE(4) - NODE(3)) * 0.5D0
         QV   = XL * (FCT(NODE(3)) + FCT(NODE(4)))
      ENDIF
C
C     store the result and return
C
      NODE(2) = QV
      IF (IERR .GT. 0) THEN
         IERR2 = 100
      ELSE
         IERR2 = 0
      ENDIF
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXINS (NMAX,DATA,STATUS,N,T,N1,N2,ITNODE,IPOS)
C[IX{GXINS}*)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE inserts two new sub-intervals into a linearly *
C  ordered list and and then re-sorts the list.                  *
C  For this a binary search is used.                             *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  NMAX  : maximum number of intervals that can be produced      *
C  DATA  : storage for the generated intervals                   *
C  N     : current number of the intervals produced              *
C  T     : pointer vector, used as an indirect address list      *
C  ITNODE: contains the status information for N1 and N2         *
C  N1, N2: sub-intervals to be inserted into DATA                *
C  STATUS: contains the status information ofr each interval     *
C  IPOS  : binary search starts at this index                    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION DATA (NMAX, 4), N1(4), N2(4)
      INTEGER          T(NMAX), STATUS(NMAX, 2), ITNODE(2)
C
C     binary search to find the position in the list
C     where the new sub-intervals are to be inserted
C
      IF (N .LE. 1) THEN
         N   = N + 2
         IT1 = T(N-1)
         IT2 = T(N)
      ELSE
         I = IPOS
         J = N
C
C        WHILE (I < J)
C
  100    K = INT((I+J)/2)
         IF ( N1(1) .LT. DATA(T(K),1) ) THEN
            I = K+1
         ELSE
            J = K-1
         ENDIF
         IF (I .LE. J) GOTO 100
C
C        insert the new sub-intervals
C
         N   = N + 2
         II  = N
         IT1 = T(N)
         IT2 = T(N-1)
         DO 30 L = N, I+2, -1
            T(L) = T(L-2)
 30      CONTINUE
C
         T(I)   = IT1
         T(I+1) = IT2
      ENDIF
C
      DO 10 J = 1, 4
         DATA(IT1, J) = N1(J)
         DATA(IT2, J) = N2(J)
 10   CONTINUE
      DO 20 J = 1, 2
         STATUS(IT1, J) = ITNODE(J)
         STATUS(IT2, J) = ITNODE(J)
 20   CONTINUE
C
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXDEL (NMAX, NCELLS, T, POS)
C[IX{GXDEL}*)
C
C*****************************************************************
C                                                                *
C  SUBROUTINE GXDEL erases the sub-interval indexed POS in the   *
C  linearly ordered list XNODES using the pointer vector T and   *
C  then re-sorts the vector T.                                   *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  NMAX  : maximum number of sub-intervals allowed               *
C  NCELLS: current number of sub-intervals produced              *
C  T     : pointer vector, used for indirect addressing          *
C  POS   : index of the interval to be erased                    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER T(NMAX), POS
C
      ISAVE  = T(POS)
      NCELLS = NCELLS - 1
C
      IF ( NCELLS .EQ. 0 ) R E T U R N
C
      DO 1 J = POS, NCELLS
         T(J) = T(J+1)
 1    CONTINUE
      T(NCELLS + 1) = ISAVE
C
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXDIV (NODE, R, L, FCT, XZERO, ZWGH, LQM, GRDFLG,
C[IX{GXDIV}*)
     +                   AK, IDGR, IERR2)
C
C*****************************************************************
C                                                                *
C  SUBROUTINE GXDIV subdivides the interval NODE into two equal  *
C  sized sub-intervals. For each of these two new sub-intervals  *
C  the numerical value of the integral, the error and the new    *
C  interval endpoints are determined.                            *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  NODE  : contains the data of the interval to be subdivided    *
C  FCT   : name of the function to be integrated                 *
C  XZERO : nodes of the GAUSS-quadratureformula                  *
C  ZWGH  : weights for the nodes                                 *
C  LQM   : quadrature formula to be used                         *
C  GRDFLG: =1, if 0 is a node; =0, if 0 is not a node            *
C  AK    : weights and nodes for the CLENSHAW-CURTIS formula     *
C  IDGR  : degree of the quadrature formula                      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  R : contains the data for the right sub-interval of NODE      *
C  L : contains the data for the left sub-interval of NODE       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : GXQUAD                                 *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      EXTERNAL         FCT
      INTEGER          IDGR
      DOUBLE PRECISION NODE(4), R(4), L(4), XZERO(IDGR),
     +                 ZWGH(IDGR), AK(0:IDGR, 2)
      INTEGER          GRDFLG
C
C     determine the new interval endpoints
C
      L(3) = NODE(3)
      R(4) = NODE(4)
      L(4) = (L(3) + R(4)) * 0.5D0
      R(3) = L(4)
C
C     integrate over the new sub-intervals
C
      CALL GXQUAD (L, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR, IERR2)
      CALL GXQUAD (R, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR, IERR2)
C
C     estimate of the error for the computed integral values
C
      L(1) = DABS(L(2) + R(2) - NODE(2))
      R(1) = L(1)
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXACC (NMAX, DATA, STATUS, N, T ,NODE, ITNODE, K)
C[IX{GXACC}*)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE permits access to the K-th sub-interval in    *
C  the linearly ordered list DATA.                               *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  NMAX  : maximum number of sub-intervals allowed               *
C  DATA  : storage for the sub-intervals produced                *
C  STATUS: contains the status information for each sub-interval *
C  N     : current number of produced sub-intervals              *
C  T     : pointer vector, used for indirect addressing          *
C  K     : determines the K-th sub-interval within DATA          *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  NODE  : contains the data for the K-th sub-interval           *
C  ITNODE: contains the status data of the K-th sub-interval     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.19.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION DATA (NMAX, 4), NODE(4)
      INTEGER          T(NMAX), STATUS (NMAX, 2), ITNODE(2)
C
      IF (K .LT. 1 .OR. K .GT. N .OR. N .GT. NMAX) RETURN
      DO 10 J = 1, 4
         NODE(J) = DATA(T(K), J)
 10   CONTINUE
      DO 20 J = 1, 2
         ITNODE(J) = STATUS(T(K), J)
 20   CONTINUE
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXGAUS (VAL, XZERO, ZWGH, IDGR, GRDFLG, FCT, QV)
C[IX{GXGAUS}*)
C
C*****************************************************************
C                                                                *
C     This SUBROUTINE determines the integral of a function FCT  *
C     over the interval [VAL(1),VAL(2)] by the GAUSS-quadrature  *
C     formula of degree IDGR                                     *
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     VAL       : vector VAL(1:2), the interval endpoints        *
C     XZERO     : vector of nodes                                *
C     ZWGH      : vector of weights                              *
C     IDGR      : number of weights and nodes                    *
C     FCT       : function to be integrated                      *
C     GRDFLG    : =1, if 0 is a node; =0, if 0 is not a node     *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     QV        : value for the integral                         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required : none                                *
C                                                                *
C*****************************************************************
C                                                                *
C     author    : Hermann-Josef Rheinbach                        *
C     editor    : Norbert Vogt                                   *
C     date      : 05.17.1989                                     *
C     source    : FORTRAN-77                                     *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER              IDGR, KDIV2, GRDFLG
      DOUBLE PRECISION     QV, XM, XL, W
      DOUBLE PRECISION     VAL(1:2), XZERO(IDGR), ZWGH(IDGR)
      EXTERNAL             FCT
C
C     GAUSS-quadrature
C
      QV    = 0.0D0
      TW    = 0.0D0
      XM    = (VAL(2) + VAL(1)) * 0.5D0
      XL    = XM - VAL(1)
      KDIV2 = INT(DBLE(IDGR) * 0.5D0)
      DO 10 I = 1, KDIV2
         W  = XL * XZERO(I)
         QV = QV + (FCT(XM-W) + FCT(XM+W)) * ZWGH(I)
 10   CONTINUE
C
C     Is 0 a node ?
C
      IF (GRDFLG .EQ. 1) THEN
         QV = (QV + FCT(XM) * ZWGH(KDIV2+1) ) * XL
      ELSE
         QV = QV * XL
      ENDIF
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXRIEP (QVAL, EL, ELEPS, EPSEXE, EPSL, M,
C[IX{GXRIEP}*)
     +                   EPS, EXE, QEXE)
C
C*****************************************************************
C                                                                *
C     SUBROUTINE GXRIEP extrapolates the current quadrature      *
C     value using the new value of T(K,0) using RICHARDSON extra-*
C     polation and the EPSILON-algorithm. This results in a new  *
C     error and a new quadrature value for the sub-interval.     *
C[BE*)
C                                                                *
C     IN/OUTPUT PARAMETERS:                                      *
C     =====================                                      *
C                                                                *
C     QVAL      : current quadrature value                       *
C     EPS       : relative accuracy                              *
C     EXE       : absolute accuracy                              *
C     EL        : vector EL(0:M), ther current ROMBERG row       *
C     ELEPS     : vector ELEPS(0:M), the current EPSILON-row     *
C     M         : number of executed calls of GXRIEP + 1         *
C     EPSEXE    : auxiliary variable                             *
C     EPSL      : auxiliary variable                             *
C     QEXE      : quadrature value                               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required : none                                *
C                                                                *
C*****************************************************************
C                                                                *
C     author    : Hermann-Josef Rheinbach                        *
C     editor    : Norbert Vogt                                   *
C     date      : 05.16.1989                                     *
C     source    : FORTRAN-77                                     *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (E)
      DOUBLE PRECISION EL(0:M), ELEPS(0:M), EPS, EPSEXE,
     +                 EPSL, EXE
      DOUBLE PRECISION SAVE1, SAVE2, AM, QEXE, QVAL
      INTEGER          M, J, IM
C
C     RICHARDSON extrapolation with the new T(K,0)
C
      IF (M .GT. 102) RETURN
      SAVE1 = EL(M)
      M     = M + 1
      AM    = 1.0D0
      EL(M) = 0.0D0
      EL1   = EL(1)
      EL(1) = QVAL
      DO 120 J = 2, M
         AM    = AM * 4.0D0
         EL2   = EL(J)
         EL(J) = (AM * EL(J-1) - EL1) / (AM-1)
         EL1   = EL2
 120  CONTINUE
C
C     extrapolation using the EPSILON algorithm
C
      IF (MOD(M,2) .EQ. 0) THEN
         IM    = M-1
         SAVE2 = ELEPS(IM)
      ELSE
         SAVE2 = ELEPS(M-2)
         IM    =  M
      ENDIF
C
      ELEPS(M) = 0.0D0
      EPSELS   = ELEPS(1)
      EPSEL1   = ELEPS(1)
      ELEPS(1) = QVAL
      DO 121 J = 2, M
         EPSEL2 = ELEPS(J)
         IF ( DABS(ELEPS(J-1)-EPSEL1) .EQ. 0.0D0) THEN
            ELEPS(J) = EPSELS
         ELSE
            ELEPS(J) = EPSELS + 1.0D0 / (ELEPS(J-1)-EPSEL1 )
            EPSELS   = EPSEL1
         ENDIF
         EPSEL1 = EPSEL2
 121  CONTINUE
      EXE = DABS(EL(M) - EL(M-1)) + DABS(EL(M) - SAVE1)
      IF (M .GT. 2) THEN
         EPSEXE = DABS(ELEPS(IM)-SAVE2) + DABS(ELEPS(IM) - ELEPS(IM-2))
         QEXE   = ELEPS(IM)
      ENDIF
      IF (EXE .LT. EPSEXE) QEXE = EL(M)
      EPSL = EPS * DABS(QEXE)
      EXE  = DMIN1 (EXE,EPSEXE)
      R E T U R N
      END

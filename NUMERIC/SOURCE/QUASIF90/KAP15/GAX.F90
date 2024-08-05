![           {Adaptive Quadrature Methods}*)                            
      SUBROUTINE GAX (INTVAL, EPS, N, FCT, NMAX, LQM, XNODES, QVAL, EXE,&
      TPOINT, STATUS, IERR, AK, IDGR)                                   
!                                                                       
!*****************************************************************      
!                                                                *      
!            Adaptive quadrature in one dimension.               *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  INTVAL: 2-dimensional array INTVAL(1:100,1:2); user supplied  *      
!          sub-intervals of integration:                         *      
!                        INTVAL(1,1) = A1, INTVAL(1,2) = B1      *      
!                           ...              ...                 *      
!                        INTVAL(N,1) = AN, INTVAL(N,2) = BN.     *      
!          Here AI must be less than BI for each I.              *      
!                                                                *      
!  EPS   : relative accuracy required for the solution           *      
!                                                                *      
!  N     : number of initially supplied sub-intervals            *      
!                                                                *      
!  FCT   : name of the function FCT(X) that is to be integrated, *      
!          provided by the user in the form                      *      
!               DOUBLE PRECISION FUNCTION FCT (X).               *      
!          In the calling program it must be defined as EXTERNAL.*      
!                                                                *      
!  NMAX  : maximal number of sub-intervals allowed to be produced*      
!          during the calculations                               *      
!                                                                *      
!  LQM   : parameter that determines which quadrature formula    *      
!          to use on the sub-intervals                           *      
!             LQM = 1: trapezoidal rule                          *      
!                 = 2: GAUSS quadrature formula                  *      
!                 = 3: CLENSHAW-CURTIS formula                   *      
!                 = 4: ROMBERG method                            *      
!                 = 5: summed NEWTON-COTES formulas              *      
!                                                                *      
!  IDGR  : this label determines the degree of the quadrature    *      
!          formula used                                          *      
!          LQM =1:IDGR not used                                  *      
!          LQM =2:IDGR = 2: GAUSS quadrature formula of degree 2 *      
!                 IDGR = 3: GAUSS quadrature formula of degree 3 *      
!                 ...                                            *      
!                 IDGR =20: GAUSS quadrature formula of degree 20*      
!                 ( 2 <= IDGR <= 20 )                            *      
!          LQM =3:IDGR = 2: CLENSHAW-CURTIS formula with         *      
!                           2 + 1 weights                        *      
!                 IDGR = 4: CLENSHAW-CURTIS formula with         *      
!                           4 + 1 weights                        *      
!                 ...                                            *      
!                 IDGR must be greater than 1 and even.          *      
!          LQM =4:IDGR not used                                  *      
!          LQM =5:IDGR = 2: summed NEWTON-COTES formula for 2    *      
!                           nodes, i.e., SIMPSON's rule          *      
!                 IDGR = 3: summed NEWTON-COTES formula for 3    *      
!                           nodes, i.e., the 3/8 formula         *      
!                 ...                                            *      
!                 IDGR = 7: summed NEWTON-COTES formula for 7    *      
!                           nodes, i.e., the 7/17280 formula     *      
!                 (We must have 2 <= IDGR <= 7 here)             *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  XNODES: vector XNODES(1:NMAX * 4) containing the sub-intervals*      
!          produced during the calculations in the following     *      
!          form:                                                 *      
!          XNODES(1) : error estimate for the quadrature from A  *      
!                      to B                                      *      
!          XNODES(2) : computed value of the integral from A to B*      
!          XNODES(3) : left hand endpoint A of the interval      *      
!          XNODES(4) : right hand endpoint B of the interval     *      
!                                                                *      
!  QVAL  : approximate value for the integral when the procedure *      
!          is stopped.                                           *      
!                                                                *      
!  EXE   : accuracy of QVAL achieved when stopping               *      
!          (absolute error)                                      *      
!                                                                *      
!  NMAX  : actual number of sub-intervals used                   *      
!                                                                *      
!  IERR  : error parameter:                                      *      
!          = 0: everything is o.k.,  EXE < EPS                   *      
!          = 1: exceeding  MAXLEV = 100 while                    *      
!               the required accuracy could not be achieved      *      
!               => choose a larger value for IDGR or a smaller   *      
!                  accuracy bound EPS                            *      
!          = 2: N > NMAX, the required accuracy was not achieved *      
!          = 3: left endpoint A > right endpoint B:              *      
!               endpoints entered incorrectly have been exchanged*      
!          = 4: conditions IERR = 3 and = 1 both hold            *      
!          = 5: conditions IERR = 3 and = 2 both hold            *      
!          = 6: degree IDGR of the GAUSS quadrature formula      *      
!               satisfies 2 > IDGR > 20                          *      
!          = 7: weights of the CLENSHAW-CURTIS formula           *      
!               IDGR are odd or less than 2                      *      
!          = 8: number of nodes IDGR in the NEWTON-COTES formula *      
!               2 > IDGR > 7                                     *      
!          = 9: N < 0  or  EPS < 1.D-12  or  NMAX < 1            *      
!          =10: 1 > LQM > 5                                      *      
!          ....                                                  *      
!                                                                *      
!          >=100: error when determining an integral value       *      
!                                                                *      
!  AUXILIARY VECTORS:                                            *      
!  ==================                                            *      
!  ITNODE: vector ITNODE(1:2); contains information on the       *      
!          current sub-interval:                                 *      
!          ITNODE(1) indicates the state of the sub-interval:    *      
!                    = 1: state = small                          *      
!                    = 0: state = big                            *      
!          ITNODE(2) contains the current level of the sub-      *      
!                    interval                                    *      
!  EL    : vector EL(1:102); contains the current row of the     *      
!          ROMBERG scheme for RICHARDSON extrapolation           *      
!  ELEPS : vector ELEPS(1:102); contains the current row of the  *      
!          EPSILON-algorithm                                     *      
!  TNODE : vector TNODE(1:4); contains the sub-interval currently*      
!          in use                                                *      
!  XZERO : vector XZERO(1:21); contains the zeros of the LEGENDRE*      
!          polynomials in ascending order                        *      
!  ZWGH  : vector ZWGH(1:21); contains the weights corresponding *      
!          to the zeros of the GAUSS-quadrature formulas         *      
!                                                                *      
!  AK    : 2-dimensional array AK(0:IDGR, 2); contains the       *      
!          weights and nodes of the CLENSHAW-CURTIS quadrature   *      
!          formula for the reference interval [-1,1] that were   *      
!          produced during the calculations                      *      
!  TPOINT: vector TPOINT(1:NMAX); used as a marker:              *      
!          this vector contains XNODES in quasi order; TPOINT(1) *      
!          points to the sub-interval with the largest quadrature*      
!          error, etc.                                           *      
!  STATUS: vector STATUS(1:NMAX * 2); contains the state of each *      
!          sub-interval (see vector ITNODE(1:4) )                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required :                                        *      
!       GXENT   : tests the input parameters                     *      
!       GXDIV   : splits the intervals                           *      
!       GXQUAD  : determines the integral approximately:         *      
!          GXGAUS  : integration using GAUSS quadrature          *      
!          CLENSH  : integration using CLENSHAW-CURTIS formulas  *      
!          QUAROM  : integration using ROMBERG method            *      
!          QUANEC  : integration using NEWTON-COTES formulas     *      
!             FCT  : function to be integrated                   *      
!       GXINS   : inserts new sub-intervals                      *      
!       GXDEL   : erases a sub-interval                          *      
!       GXACC   : accesses a sub-interval                        *      
!       GXRIEP  : RICHARDSON-extrapolation, EPSILON-algorithm    *      
!       GALE0   : nodes and weights for GAUSS-quadrature         *      
!          GXPOLY  : HORNER scheme                               *      
!          GXPEGA  : PEGASUS method                              *      
!       WGKNOT  : weights, nodes for CLENSHAW-CURTIS formulas    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  sources : 1. Eul, T. and  Rheinbach, H.J.:                    *      
!               Lecture: The combination of adaptive and extra-  *      
!               polation methods for numerical integration,      *      
!               RWTH Aachen, October 1984.                       *      
!            2. D. Kahaner and J. Stoer, see [KAHA83].           *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  editor   : Norbert Vogt                                       *      
!  date     : 04.10.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      EXTERNAL FCT 
      INTEGER IDGR, GRDFLG, IERR2 
      DOUBLEPRECISION TNODE (4), N1 (4), N2 (4), XNODES (1:NMAX, 1:4) 
      DOUBLEPRECISION EL (0:102), ELEPS (0:102), INTVAL (1:100, 1:2) 
      DOUBLEPRECISION XZERO (21), ZWGH (21), AK (0:IDGR, 2), HELP, EPS 
      INTEGER TPOINT (NMAX), STATUS (1:NMAX, 1:2), ITNODE (2) 
      DOUBLEPRECISION QVAL, EXE, QEXE, EPSL, EPSEXE, BEE, TEE 
!                                                                       
!     testing the input parameters,                                     
!     if necessary return                                               
!                                                                       
      IERR = 0 
      CALL GXENT (LQM, N, NMAX, EPS, AK, IDGR, ZWGH, XZERO, GRDFLG,     &
      IERR)                                                             
      IF (IERR.GT.0) RETURN 
!                                                                       
!     initialize the pointer vector TPOINT                              
!                                                                       
      DO 10 I = 1, NMAX 
         TPOINT (I) = I 
   10 END DO 
!                                                                       
!     initialize the vector XNODES for future calculations              
!                                                                       
      QVAL = 0.0D0 
      TEE = 0.0D0 
      TNODE (1) = 1000.0D0 
      DO 30 I = 1, N 
!                                                                       
!        the integral is computed for the I-th                          
!        starting interval                                              
!                                                                       
         TNODE (3) = INTVAL (I, 1) 
         TNODE (4) = INTVAL (I, 2) 
!                                                                       
!        check that left endpoint < right endpoint                      
!                                                                       
         IF (DABS (TNODE (3) ) .GT.DABS (TNODE (4) ) ) THEN 
!                                                                       
!           if right endpoint < left endpoint:                          
!           ==> exchange the two                                        
!                                                                       
            IERR = 3 
            HELP = TNODE (3) 
            TNODE (3) = TNODE (4) 
            TNODE (4) = HELP 
         ENDIF 
!                                                                       
         CALL GXQUAD (TNODE, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR,   &
         IERR2)                                                         
!                                                                       
!        store the I-th sub-interval in the vector XNODES               
!                                                                       
         DO 40 J = 1, 4 
            XNODES (I, J) = TNODE (J) 
   40    END DO 
!                                                                       
         STATUS (I, 1) = 1 
         STATUS (I, 2) = 0 
         QVAL = QVAL + TNODE (2) 
         TEE = TEE+1000.0D0 
   30 END DO 
!                                                                       
!     initialize remaining variables                                    
!                                                                       
      LEVEL = 0 
      MAXLEV = 100 
      BEE = 0.0D0 
      EXE = TEE 
      EPSEXE = TEE 
      EL (1) = QVAL 
      ELEPS (1) = QVAL 
      EPSL = 0.0D0 
      QEXE = 1000.0D0 
      M = 1 
!                                                                       
!     integrate                                                         
!                                                                       
      DO 111 ITER = 1, NMAX 
!                                                                       
!        test break-off criterion of the method used                    
!                                                                       
         IF (TEE.LE.EPSL.OR.EXE.LE.EPSL.OR.LEVEL.GT.MAXLEV) THEN 
!                                                                       
!           If TEE <= EPSL  or  EXE <= EPSL  or  LEVEL > MAXLEV         
!           ==> correct computation                                     
!                                                                       
            IF (IERR2.GT.0) IERR = IERR + IERR2 
            IF (LEVEL.GT.MAXLEV) IERR = IERR + 1 
            IF (EXE.LT.TEE) THEN 
               QVAL = QEXE 
            ELSEIF (TEE.LT.EXE) THEN 
               EXE = TEE 
            ENDIF 
            NMAX = N 
            RETURN 
         ENDIF 
!                                                                       
!        retrieve the interval with the largest error and transfer      
!        it to the vector TNODE                                         
!                                                                       
         CALL GXACC (NMAX, XNODES, STATUS, N, TPOINT, TNODE, ITNODE, 1) 
         CALL GXDEL (NMAX, N, TPOINT, 1) 
         CALL GXDIV (TNODE, N1, N2, FCT, XZERO, ZWGH, LQM, GRDFLG, AK,  &
         IDGR, IERR2)                                                   
!                                                                       
!                                                                       
         IF (ITNODE (2) .EQ.LEVEL) THEN 
            LEVEL = LEVEL + 1 
!                                                                       
!           all used sub-intervals are assigned the state 'BIG'.        
!           This requires to recalculate the value for BEE              
!                                                                       
            BEE = 0.0D0 
            DO 100 I = 1, N 
               STATUS (TPOINT (I), 1) = 0 
               BEE = BEE+XNODES (TPOINT (I), 1) 
  100       END DO 
            ITNODE (1) = 1 
            ITNODE (2) = LEVEL 
         ELSE 
            IF (ITNODE (2) .EQ. (LEVEL - 1) ) THEN 
               ITNODE (1) = 1 
               BEE = DABS (BEE-TNODE (1) ) 
            ELSE 
               ITNODE (1) = 0 
               BEE = DABS (BEE+ (N1 (1) + N2 (1) - TNODE (1) ) ) 
            ENDIF 
            ITNODE (2) = ITNODE (2) + 1 
         ENDIF 
         IPOS = 1 
!                                                                       
!        check whether two additional sub-intervals can be inserted     
!                                                                       
    2    IF ( (N + 2) .GT.NMAX) THEN 
!                                                                       
!           If N+2 > NMAX                                               
!           ==>  correct computations                                   
!                                                                       
            IERR = IERR + 2 
            IF (IERR2.GT.0) IERR = IERR + IERR2 
            IF (EXE.LT.TEE) THEN 
               QVAL = QEXE 
            ELSEIF (TEE.LT.EXE) THEN 
               EXE = TEE 
            ENDIF 
            NMAX = N 
            RETURN 
         ENDIF 
         CALL GXINS (NMAX, XNODES, STATUS, N, TPOINT, N1, N2, ITNODE,   &
         IPOS)                                                          
!                                                                       
!        updating of values  for TEE  and  QVAL                         
!                                                                       
         TEE = DABS (TEE+ (N1 (1) + N2 (1) - TNODE (1) ) ) 
         QVAL = QVAL + (N1 (2) + N2 (2) - TNODE (2) ) 
!                                                                       
         IF (BEE.GT.EPSL.AND.TEE.GT.EPSL) THEN 
!                                                                       
!           find the sub-interval with state 'BIG'                      
!           that has the largest quadrature error                       
!                                                                       
            DO 110 I = IPOS, N 
               IF (STATUS (TPOINT (I), 1) .NE.1) THEN 
                  IPOS = I 
                  CALL GXACC (NMAX, XNODES, STATUS, N, TPOINT, TNODE,   &
                  ITNODE, IPOS)                                         
                  CALL GXDEL (NMAX, N, TPOINT, IPOS) 
                  CALL GXDIV (TNODE, N1, N2, FCT, XZERO, ZWGH, LQM,     &
                  GRDFLG, AK, IDGR, IERR2)                              
!                                                                       
!                 are the new rows of the state 'BIG' or 'SMALL' ?      
!                                                                       
                  IF (ITNODE (2) .EQ. (LEVEL - 1) ) THEN 
                     ITNODE (1) = 1 
                     BEE = DABS (BEE-TNODE (1) ) 
                  ELSE 
                     ITNODE (1) = 0 
                     BEE = DABS (BEE+ (N1 (1) + N2 (1) - TNODE (1) ) ) 
                  ENDIF 
                  ITNODE (2) = ITNODE (2) + 1 
!                                                                       
                  GOTO 2 
!                                                                       
               ENDIF 
  110       END DO 
         ENDIF 
         CALL GXRIEP (QVAL, EL, ELEPS, EPSEXE, EPSL, M, EPS, EXE, QEXE) 
  111 END DO 
      IF (IERR2.GT.0) IERR = IERR + IERR2 
      IF (EXE.LT.TEE) THEN 
         QVAL = QEXE 
      ELSEIF (TEE.LT.EXE) THEN 
         EXE = TEE 
      ENDIF 
      NMAX = ITER 
      RETURN 
      END SUBROUTINE GAX                            
!                                                                       
!                                                                       
      SUBROUTINE GXENT (LQM, N, NMAX, EPS, AK, IDGR, ZWGH, XZERO,       &
      GRDFLG, IERR)                                                     
!                                                                       
!*****************************************************************      
!                                                                *      
!     This SUBROUTINE tests the validity of the input parameters,*      
!     and, if necessary, determines the weights for a CLENSHAW-  *      
!     CURTIS or GAUSS quadrature formula of degree IDGR.         *      
!                                                                *      
!                                                                *      
!     INPUT PARAMERTERS:                                         *      
!     ==================                                         *      
!     LQM       : quadrature formula to be used                  *      
!     N         : number of starting sub-intervals               *      
!     NMAX      : maximum number of sub-intervals that can be    *      
!                 produced                                       *      
!     EPS       : relative accuracy bound                        *      
!     IDGR      : degree of the quadrature formula               *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     AK        : 2-dimensional array AK(0:IDGR, 2), the weights *      
!                 and nodes for CLENSHAW-CURTIS                  *      
!     ALPHA     : vector ALPHA(1:21), the nodes for GAUSS        *      
!     ZWGH      : vector ZWGH(1:21), the weights for GAUSS       *      
!     XZERO     : vector XZERO(1:21), the zeros for GAUSS        *      
!     GRDFLG    : =1, if 0 is a node; =0, if zero is not a node  *      
!     IERR      : error parameter                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required : GALE0, WGKNOT                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     author   : Norbert Vogt                                    *      
!     date     : 05.18.1989                                      *      
!     source   : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, IDGR, IERR, GRDFLG, NMAX 
      DOUBLEPRECISION AK (0:IDGR, 2), ZWGH (21), XZERO (21) 
      DOUBLEPRECISION EPS 
!                                                                       
!     testing the input parameter                                       
!                                                                       
      IF (N.LT.0.OR.EPS.LT.1.0D-12.OR.NMAX.LT.1) THEN 
!                                                                       
!        If  N < 0   or  EPS < 1.D-12  or  NMAX < 1                     
!        ==> return to calling program                                  
!                                                                       
         IERR = 9 
         RETURN 
      ELSEIF (LQM.LT.1.OR.LQM.GT.5) THEN 
!                                                                       
!        If  LQM < 1  or  LQM > 5                                       
!        ==> return to calling program                                  
!                                                                       
         IERR = 10 
         RETURN 
      ELSEIF (LQM.EQ.2) THEN 
!                                                                       
!        If LQM = 2                                                     
!        ==> use GAUSS quadrature                                       
!                                                                       
         IF (IDGR.LE.1.OR.IDGR.GT.20) THEN 
!                                                                       
!           If  IDGR <= 1  or  IDGR > 20                                
!           ==> return to calling program                               
!                                                                       
            IERR = 6 
            RETURN 
         ELSE 
!                                                                       
!           the desired GAUSS quadrature formula                        
!           of degree 1 < IDGR  < 20 is generated                       
!                                                                       
            CALL GALE0 (IDGR, .TRUE., XZERO, ZWGH) 
!                                                                       
!           nodes                                                       
!                                                                       
            GRDFLG = 0 
!                                                                       
!           Is 0 a node ?                                               
!                                                                       
            IF ( (MOD (IDGR, 2) .EQ.1) ) GRDFLG = 1 
         ENDIF 
      ELSEIF (LQM.EQ.3) THEN 
!                                                                       
!        If  LQM = 3                                                    
!        ==> use CLENSHAW-CURTIS formulas                               
!                                                                       
         IF (MOD (IDGR, 2) .NE.0.OR.IDGR.LT.2) THEN 
!                                                                       
!           If  IDGR is not even or if IDGR < 2                         
!           ==> return to calling program                               
!                                                                       
            IERR = 7 
            RETURN 
         ELSE 
!                                                                       
!           the weights and nodes for the CLENSHAW-                     
!           CURTIS formulas are determined                              
!                                                                       
            CALL WGKNOT (IDGR, AK, IERR2) 
         ENDIF 
      ELSEIF (LQM.EQ.5) THEN 
!                                                                       
!        If  LQM = 5                                                    
!        ==> use summed NEWTON-COTES formula                            
!                                                                       
         IF (IDGR.LT.2.OR.IDGR.GT.7) THEN 
!                                                                       
!           If  IDGR < 2 or IDGR > 7                                    
!           ==> return to calling program                               
!                                                                       
            IERR = 8 
            RETURN 
         ENDIF 
!                                                                       
!        for LQM = 1 and LQM = 4 no further                             
!        tests are required                                             
!                                                                       
      ENDIF 
      RETURN 
      END SUBROUTINE GXENT                          
!                                                                       
!                                                                       
      SUBROUTINE GXQUAD (NODE, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR, &
      IERR2)                                                            
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE determines the approximate integral of the    *      
!  FUNCTION FCT(X) over the interval (A,B) by the quadrature     *      
!  formula determined by LQM as follows:                         *      
!        LQM = 1:  trapezoidal rule                              *      
!        LQM = 2:  IDGR =2,...,20 : I-point-formula of GAUSS     *      
!        LQM = 3:  IDGR =2,4,6,...: CLENSHAW-CURTIS formula with *      
!                                   IDGR + 1 weights             *      
!        LQM = 4:  ROMBERG method                                *      
!        LQM = 5:  IDGR =2,...,7  : NEWTON-COTES formulas with   *      
!                                   IDGR nodes                   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NODE  : contains the index of the sub-interval over which we  *      
!          integrate                                             *      
!  FCT   : name of integrand                                     *      
!  XZERO : nodes of the GAUSS-quadrature formula                 *      
!  ZWGH  : weights for the nodes of the GAUSS-formulas           *      
!  LQM   : indicates quadrature formula to be used               *      
!  GRDFLG: =1, if 0 is a node; =0, if 0 is not a node            *      
!  AK    : weights and nodes for the CLENSHAW-CURTIS formulas    *      
!  IDGR  : degree of the quadrature formula                      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NODE(2): computed value for the integral over the interval    *      
!           [NODE(3), NODE(4)]                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : CLENSH, QUAROM, QUANEC, GXGAUSS        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  editor   : Norbert Vogt                                       *      
!  date     : 05.10.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      EXTERNAL FCT 
      INTEGER IDGR 
      DOUBLEPRECISION NODE (4), XZERO (IDGR), ZWGH (IDGR), XL 
      DOUBLEPRECISION AK (0:IDGR, 2), H, EL (10), ERREST, QV 
      INTEGER GRDFLG 
!                                                                       
      IF (LQM.EQ.2) THEN 
!                                                                       
!        If  LQM = 2                                                    
!        ==> the integral is determined using the GAUSS-quadrature      
!            formula of degree IDGR                                     
!                                                                       
         CALL GXGAUS (NODE (3), XZERO, ZWGH, IDGR, GRDFLG, FCT, QV) 
!                                                                       
      ELSEIF (LQM.EQ.3) THEN 
!                                                                       
!        If  LQM = 2                                                    
!        ==> the integral is determined using the summed CLENSHAW-      
!            CURTIS formula with IDGR weights                           
!                                                                       
         CALL CLENSH (FCT, IDGR, NODE (3), 1, AK, QV, IERR) 
!                                                                       
      ELSEIF (LQM.EQ.4) THEN 
!                                                                       
!        If LQM = 4                                                     
!        ==> the integral is determined via the ROMBERG method          
!                                                                       
         H = 0.0D0 
         NROWS = 10 
         CALL QUAROM (NODE (3), NODE (4), 1.0D-8, NROWS, H, FCT, EL, QV,&
         ERREST, IERR)                                                  
!                                                                       
      ELSEIF (LQM.EQ.5) THEN 
!                                                                       
!        If  LQM = 5                                                    
!        ==> the integral is determined by the summed                   
!            NEWTON-COTES formula with IDGR nodes                       
!                                                                       
         CALL QUANEC (NODE (3), NODE (4), 1, IDGR, FCT, QV, EL (1),     &
         EL (2), IERR)                                                  
!                                                                       
      ELSE 
!                                                                       
!        If  LQM = 1                                                    
!        ==> the integral is determined by the trapezoidal rule         
!                                                                       
         XL = (NODE (4) - NODE (3) ) * 0.5D0 
         QV = XL * (FCT (NODE (3) ) + FCT (NODE (4) ) ) 
      ENDIF 
!                                                                       
!     store the result and return                                       
!                                                                       
      NODE (2) = QV 
      IF (IERR.GT.0) THEN 
         IERR2 = 100 
      ELSE 
         IERR2 = 0 
      ENDIF 
      RETURN 
      END SUBROUTINE GXQUAD                         
!                                                                       
!                                                                       
      SUBROUTINE GXINS (NMAX, DATA, STATUS, N, T, N1, N2, ITNODE, IPOS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE inserts two new sub-intervals into a linearly *      
!  ordered list and and then re-sorts the list.                  *      
!  For this a binary search is used.                             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NMAX  : maximum number of intervals that can be produced      *      
!  DATA  : storage for the generated intervals                   *      
!  N     : current number of the intervals produced              *      
!  T     : pointer vector, used as an indirect address list      *      
!  ITNODE: contains the status information for N1 and N2         *      
!  N1, N2: sub-intervals to be inserted into DATA                *      
!  STATUS: contains the status information ofr each interval     *      
!  IPOS  : binary search starts at this index                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  editor   : Norbert Vogt                                       *      
!  date     : 04.10.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION DATA (NMAX, 4), N1 (4), N2 (4) 
      INTEGER T (NMAX), STATUS (NMAX, 2), ITNODE (2) 
!                                                                       
!     binary search to find the position in the list                    
!     where the new sub-intervals are to be inserted                    
!                                                                       
      IF (N.LE.1) THEN 
         N = N + 2 
         IT1 = T (N - 1) 
         IT2 = T (N) 
      ELSE 
         I = IPOS 
         J = N 
!                                                                       
!        WHILE (I < J)                                                  
!                                                                       
  100    K = INT ( (I + J) / 2) 
         IF (N1 (1) .LT.DATA (T (K), 1) ) THEN 
            I = K + 1 
         ELSE 
            J = K - 1 
         ENDIF 
         IF (I.LE.J) GOTO 100 
!                                                                       
!        insert the new sub-intervals                                   
!                                                                       
         N = N + 2 
         II = N 
         IT1 = T (N) 
         IT2 = T (N - 1) 
         DO 30 L = N, I + 2, - 1 
            T (L) = T (L - 2) 
   30    END DO 
!                                                                       
         T (I) = IT1 
         T (I + 1) = IT2 
      ENDIF 
!                                                                       
      DO 10 J = 1, 4 
         DATA (IT1, J) = N1 (J) 
         DATA (IT2, J) = N2 (J) 
   10 END DO 
      DO 20 J = 1, 2 
         STATUS (IT1, J) = ITNODE (J) 
         STATUS (IT2, J) = ITNODE (J) 
   20 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE GXINS                          
!                                                                       
!                                                                       
      SUBROUTINE GXDEL (NMAX, NCELLS, T, POS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  SUBROUTINE GXDEL erases the sub-interval indexed POS in the   *      
!  linearly ordered list XNODES using the pointer vector T and   *      
!  then re-sorts the vector T.                                   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NMAX  : maximum number of sub-intervals allowed               *      
!  NCELLS: current number of sub-intervals produced              *      
!  T     : pointer vector, used for indirect addressing          *      
!  POS   : index of the interval to be erased                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  editor   : Norbert Vogt                                       *      
!  date     : 04.10.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER T (NMAX), POS 
!                                                                       
      ISAVE = T (POS) 
      NCELLS = NCELLS - 1 
!                                                                       
      IF (NCELLS.EQ.0) RETURN 
!                                                                       
      DO 1 J = POS, NCELLS 
         T (J) = T (J + 1) 
    1 END DO 
      T (NCELLS + 1) = ISAVE 
!                                                                       
      RETURN 
      END SUBROUTINE GXDEL                          
!                                                                       
!                                                                       
      SUBROUTINE GXDIV (NODE, R, L, FCT, XZERO, ZWGH, LQM, GRDFLG, AK,  &
      IDGR, IERR2)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  SUBROUTINE GXDIV subdivides the interval NODE into two equal  *      
!  sized sub-intervals. For each of these two new sub-intervals  *      
!  the numerical value of the integral, the error and the new    *      
!  interval endpoints are determined.                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NODE  : contains the data of the interval to be subdivided    *      
!  FCT   : name of the function to be integrated                 *      
!  XZERO : nodes of the GAUSS-quadratureformula                  *      
!  ZWGH  : weights for the nodes                                 *      
!  LQM   : quadrature formula to be used                         *      
!  GRDFLG: =1, if 0 is a node; =0, if 0 is not a node            *      
!  AK    : weights and nodes for the CLENSHAW-CURTIS formula     *      
!  IDGR  : degree of the quadrature formula                      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  R : contains the data for the right sub-interval of NODE      *      
!  L : contains the data for the left sub-interval of NODE       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : GXQUAD                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  editor   : Norbert Vogt                                       *      
!  date     : 04.10.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      EXTERNAL FCT 
      INTEGER IDGR 
      DOUBLEPRECISION NODE (4), R (4), L (4), XZERO (IDGR), ZWGH (IDGR),&
      AK (0:IDGR, 2)                                                    
      INTEGER GRDFLG 
!                                                                       
!     determine the new interval endpoints                              
!                                                                       
      L (3) = NODE (3) 
      R (4) = NODE (4) 
      L (4) = (L (3) + R (4) ) * 0.5D0 
      R (3) = L (4) 
!                                                                       
!     integrate over the new sub-intervals                              
!                                                                       
      CALL GXQUAD (L, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR, IERR2) 
      CALL GXQUAD (R, FCT, XZERO, ZWGH, LQM, GRDFLG, AK, IDGR, IERR2) 
!                                                                       
!     estimate of the error for the computed integral values            
!                                                                       
      L (1) = DABS (L (2) + R (2) - NODE (2) ) 
      R (1) = L (1) 
      RETURN 
      END SUBROUTINE GXDIV                          
!                                                                       
!                                                                       
      SUBROUTINE GXACC (NMAX, DATA, STATUS, N, T, NODE, ITNODE, K) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE permits access to the K-th sub-interval in    *      
!  the linearly ordered list DATA.                               *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NMAX  : maximum number of sub-intervals allowed               *      
!  DATA  : storage for the sub-intervals produced                *      
!  STATUS: contains the status information for each sub-interval *      
!  N     : current number of produced sub-intervals              *      
!  T     : pointer vector, used for indirect addressing          *      
!  K     : determines the K-th sub-interval within DATA          *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NODE  : contains the data for the K-th sub-interval           *      
!  ITNODE: contains the status data of the K-th sub-interval     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  editor   : Norbert Vogt                                       *      
!  date     : 04.19.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION DATA (NMAX, 4), NODE (4) 
      INTEGER T (NMAX), STATUS (NMAX, 2), ITNODE (2) 
!                                                                       
      IF (K.LT.1.OR.K.GT.N.OR.N.GT.NMAX) RETURN 
      DO 10 J = 1, 4 
         NODE (J) = DATA (T (K), J) 
   10 END DO 
      DO 20 J = 1, 2 
         ITNODE (J) = STATUS (T (K), J) 
   20 END DO 
      RETURN 
      END SUBROUTINE GXACC                          
!                                                                       
!                                                                       
      SUBROUTINE GXGAUS (VAL, XZERO, ZWGH, IDGR, GRDFLG, FCT, QV) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This SUBROUTINE determines the integral of a function FCT  *      
!     over the interval [VAL(1),VAL(2)] by the GAUSS-quadrature  *      
!     formula of degree IDGR                                     *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     VAL       : vector VAL(1:2), the interval endpoints        *      
!     XZERO     : vector of nodes                                *      
!     ZWGH      : vector of weights                              *      
!     IDGR      : number of weights and nodes                    *      
!     FCT       : function to be integrated                      *      
!     GRDFLG    : =1, if 0 is a node; =0, if 0 is not a node     *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     QV        : value for the integral                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required : none                                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     author    : Hermann-Josef Rheinbach                        *      
!     editor    : Norbert Vogt                                   *      
!     date      : 05.17.1989                                     *      
!     source    : FORTRAN-77                                     *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER IDGR, KDIV2, GRDFLG 
      DOUBLEPRECISION QV, XM, XL, W 
      DOUBLEPRECISION VAL (1:2), XZERO (IDGR), ZWGH (IDGR) 
      EXTERNAL FCT 
!                                                                       
!     GAUSS-quadrature                                                  
!                                                                       
      QV = 0.0D0 
      TW = 0.0D0 
      XM = (VAL (2) + VAL (1) ) * 0.5D0 
      XL = XM - VAL (1) 
      KDIV2 = INT (DBLE (IDGR) * 0.5D0) 
      DO 10 I = 1, KDIV2 
         W = XL * XZERO (I) 
         QV = QV + (FCT (XM - W) + FCT (XM + W) ) * ZWGH (I) 
   10 END DO 
!                                                                       
!     Is 0 a node ?                                                     
!                                                                       
      IF (GRDFLG.EQ.1) THEN 
         QV = (QV + FCT (XM) * ZWGH (KDIV2 + 1) ) * XL 
      ELSE 
         QV = QV * XL 
      ENDIF 
      RETURN 
      END SUBROUTINE GXGAUS                         
!                                                                       
!                                                                       
      SUBROUTINE GXRIEP (QVAL, EL, ELEPS, EPSEXE, EPSL, M, EPS, EXE,    &
      QEXE)                                                             
!                                                                       
!*****************************************************************      
!                                                                *      
!     SUBROUTINE GXRIEP extrapolates the current quadrature      *      
!     value using the new value of T(K,0) using RICHARDSON extra-*      
!     polation and the EPSILON-algorithm. This results in a new  *      
!     error and a new quadrature value for the sub-interval.     *      
!                                                                *      
!     IN/OUTPUT PARAMETERS:                                      *      
!     =====================                                      *      
!                                                                *      
!     QVAL      : current quadrature value                       *      
!     EPS       : relative accuracy                              *      
!     EXE       : absolute accuracy                              *      
!     EL        : vector EL(0:M), ther current ROMBERG row       *      
!     ELEPS     : vector ELEPS(0:M), the current EPSILON-row     *      
!     M         : number of executed calls of GXRIEP + 1         *      
!     EPSEXE    : auxiliary variable                             *      
!     EPSL      : auxiliary variable                             *      
!     QEXE      : quadrature value                               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required : none                                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     author    : Hermann-Josef Rheinbach                        *      
!     editor    : Norbert Vogt                                   *      
!     date      : 05.16.1989                                     *      
!     source    : FORTRAN-77                                     *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (E) 
      DOUBLEPRECISION EL (0:M), ELEPS (0:M), EPS, EPSEXE, EPSL, EXE 
      DOUBLEPRECISION SAVE1, SAVE2, AM, QEXE, QVAL 
      INTEGER M, J, IM 
!                                                                       
!     RICHARDSON extrapolation with the new T(K,0)                      
!                                                                       
      IF (M.GT.102) RETURN 
      SAVE1 = EL (M) 
      M = M + 1 
      AM = 1.0D0 
      EL (M) = 0.0D0 
      EL1 = EL (1) 
      EL (1) = QVAL 
      DO 120 J = 2, M 
         AM = AM * 4.0D0 
         EL2 = EL (J) 
         EL (J) = (AM * EL (J - 1) - EL1) / (AM - 1) 
         EL1 = EL2 
  120 END DO 
!                                                                       
!     extrapolation using the EPSILON algorithm                         
!                                                                       
      IF (MOD (M, 2) .EQ.0) THEN 
         IM = M - 1 
         SAVE2 = ELEPS (IM) 
      ELSE 
         SAVE2 = ELEPS (M - 2) 
         IM = M 
      ENDIF 
!                                                                       
      ELEPS (M) = 0.0D0 
      EPSELS = ELEPS (1) 
      EPSEL1 = ELEPS (1) 
      ELEPS (1) = QVAL 
      DO 121 J = 2, M 
         EPSEL2 = ELEPS (J) 
         IF (DABS (ELEPS (J - 1) - EPSEL1) .EQ.0.0D0) THEN 
            ELEPS (J) = EPSELS 
         ELSE 
            ELEPS (J) = EPSELS + 1.0D0 / (ELEPS (J - 1) - EPSEL1) 
            EPSELS = EPSEL1 
         ENDIF 
         EPSEL1 = EPSEL2 
  121 END DO 
      EXE = DABS (EL (M) - EL (M - 1) ) + DABS (EL (M) - SAVE1) 
      IF (M.GT.2) THEN 
         EPSEXE = DABS (ELEPS (IM) - SAVE2) + DABS (ELEPS (IM) - ELEPS (&
         IM - 2) )                                                      
         QEXE = ELEPS (IM) 
      ENDIF 
      IF (EXE.LT.EPSEXE) QEXE = EL (M) 
      EPSL = EPS * DABS (QEXE) 
      EXE = DMIN1 (EXE, EPSEXE) 
      RETURN 
      END SUBROUTINE GXRIEP                         

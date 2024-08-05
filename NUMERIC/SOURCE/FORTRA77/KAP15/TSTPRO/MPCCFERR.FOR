      PROGRAM TEST
C****************************************************************************
C                                                                           *
C  Test program for the subroutines CCFERR, CLENSH und WGKNOT.              *
C                                                                           *
C  We compute the integral I(FCT; A,B), for  FCT=F(X)=X*LN(X), A=1, B=4,    *
C  using the CLenshaw-Curtis formula with global error order  N+2 for N=2,4 *
C  and 6.                                                                   *
C  The integrand FCT appears as a DOUBLE PRECISION FUNCTION after the       *
C  program.                                                                 *
C  Results for N=2,4,6 :                                                    *
C                                                                           *
C[                                                                         ]*
C[                                                                         ]*
C[ COMPUTING AN INTEGRAL WITH THE CLENSHAW-CURTIS FORMULA:                 ]*
C[ =======================================================                 ]*
C[ COMPUTE  I(FCT; A,B)                                                    ]*
C[ FOR: FCT = F(X) = X*LN(X)                                               ]*
C[      A = 1, B = 4                                                       ]*
C[ PARTITION USED Z: 1, 1.5, 2, 2.5, 3, 3.5, 4                             ]*
C[                                                                         ]*
C[ EXACT SOLUTION: I(FCT; A,B) =  7.34035488896                            ]*
C[                                                                         ]*
C[                                                                         ]*
C[ OUTPUT:                                                                 ]*
C[ =======                                                                 ]*
C[ N      = SIZE OF PARTITION, NODES - 1                                   ]*
C[ QCCNZ  = RESULT USING PARTITION Z                                       ]*
C[ QCCNZ2 = RESULT USING PARTITION Z/2                                     ]*
C[ DIFF1  = DIFFERENCE QCCNZ AND QCCNZ2                                    ]*
C[ SCHETZ = ESTIMATE FOR GLOBAL PROCEDURAL ERROR OF THE RESULT FOR Z/2     ]*
C[ QCCNST = IMPROVED APPROXIMATE VALUE FROM QCCNZ2                         ]*
C[ DIFF2  = DIFFERENCE BETWEEN EXACT SOLUTION AND QCCNZ2                   ]*
C[                                                                         ]*
C[                                                                         ]*
C[ N    QCCNZ       QCCNZ2     DIFF1      SCHETZ      QCCNST      DIFF2    ]*
C[ ----------------------------------------------------------------------- ]*
C[  2  7.34037437  7.34035615 .1823E-04 -.12151E-05  7.34035493  .4190E-07 ]*
C[  4  7.34035487  7.34035489 .1621E-07  .25728E-09  7.34035489  .3140E-10 ]*
C[  6  7.34035489  7.34035489 .1388E-10  .54426E-13  7.34035489  .9770E-14 ]*
C[                                                                         ]*
C[                                                                         ]*
C                                                                           *
C****************************************************************************
C
C  Declarations
C
      PARAMETER (M=6,NV=6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION Z(0:M),WORK(0:NV,2),Z2(0:2*M)
      EXTERNAL FCT
C
C  set up partition
C
      DATA (Z(I),I=0,M) /1.0D0,1.5D0,2.0D0,2.5D0,3.0D0,3.5D0,4.0D0/
C
C  Exact solution I(FCT; A,B) from theory,
C  for A=1, B=4, FCT=X*LN(X).
C
      EL = 8.D0*DLOG(4.D0)-3.75D0
      WRITE (*,1000) EL
      WRITE (*,1001)
      WRITE (*,1002)
C
C  Compute the integral
C
      N = 0
      DO 100 N = 2,6,2
      CALL CCFERR (FCT,N,Z,Z2,M,WORK,QCCNZ,QCCNZ2,SCHETZ,QCCNST,IFEHL)
      IF (IFEHL .NE. 0) THEN
         WRITE (*,*)'ERROR. IFEHL=',IFEHL,' ,N=',N,' ,M=',M
         STOP
      END IF
      DIFF1 = DABS(QCCNZ-QCCNZ2)
      DIFF2 = DABS(EL-QCCNST)
      WRITE (*,1100)N,QCCNZ,QCCNZ2,DIFF1,SCHETZ,QCCNST,DIFF2
  100 CONTINUE
      WRITE (*,1200)
      STOP
 1000 FORMAT(1X,'C[',T77,']*',/,
     +       1X,'C[',T77,']*',/,
     +       1X,'C[ COMPUTING AN INTEGRAL WITH THE CLENSHAW-CURTIS ',
     +          'FORMULA:',T77,']*',/,
     +       1X,'C[ ',55('='),T77,']*',/,
     +       1X,'C[ COMPUTE  I(FCT; A,B)',T77,']*',/,
     +       1X,'C[ FOR: FCT = F(X) = X*LN(X)',T77,']*',/,
     +       1X,'C[      A = 1, B = 4',T77,']*',/,
     +       1X,'C[ PARTITION USED Z: 1, 1.5, 2, 2.5, 3, ',
     +          '3.5, 4',T77,']*',/,
     +       1X,'C[',T77,']*',/,
     +       1X,'C[ EXACT SOLUTION: I(FCT; A,B) = ',F14.11,T77,']*',/,
     +       1X,'C[',T77,']*',/,
     +       1X,'C[',T77,']*',/,
     +       1X,'C[ OUTPUT:',T77,']*',/,
     +       1X,'C[ ',7('='),T77,']*',/,
     +       1X,'C[ N      = SIZE OF PARTITION, NODES - 1 ',T77,']*',/,
     +       1X,'C[ QCCNZ  = RESULT USING PARTITION Z ',T77,']*')
 1001 FORMAT(1X,'C[ QCCNZ2 = RESULT USING PARTITION Z/2 ',T77,']*',/,
     +       1X,'C[ DIFF1  = DIFFERENCE QCCNZ AND QCCNZ2 ',T77,']*',/,
     +       1X,'C[ SCHETZ = ESTIMATE FOR GLOBAL PROCEDURAL ERROR OF ',
     +          'THE RESULT FOR Z/2 ',T77,']*',/,
     +       1X,'C[ QCCNST = IMPROVED APPROXIMATE VALUE FROM QCCNZ2',
     +          T77,']*',/,
     +       1X,'C[ DIFF2  = DIFFERENCE BETWEEN EXACT SOLUTION AND ',
     +          'QCCNZ2',T77,']*',/,
     +       1X,'C[',T77,']*',/,
     +       1X,'C[',T77,']*')
 1002 FORMAT(1X,'C[ N',4X,'QCCNZ',7X,'QCCNZ2',5X,'DIFF1',6X,'SCHETZ',6X,
     +          'QCCNST',6X,'DIFF2',T77,']*',/,
     +       1X,'C[ ',71('-'),T77,']*')
 1100 FORMAT(1X,'C[ ',I2,2(2X,F10.8),1X,E9.4,1X,E11.5,2X,F10.8,2X,E9.4,
     +          T77,']*')
 1200 FORMAT(1X,'C[',T77,']*',/,
     +       1X,'C[',T77,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION FCT(X)
C
C***********************************************************************
C  Integrand
C***********************************************************************
C
      DOUBLE PRECISION X
      FCT = X*DLOG(X)
      RETURN
      END

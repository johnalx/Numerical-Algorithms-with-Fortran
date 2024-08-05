      PROGRAM TEST
C***************************************************************************
C                                                                          *
C  Test program for the subroutines CURVPG, KNOTVG and DEBOOR.             *
C  To compute a uniform closed B spline.                                   *
C  Berechnung einer uniformen geschlossenen B-Spline-KURVE.                *
C                                                                          *
C--------------------------------------------------------------------------*
C                                                                          *
C     The computed results are shown with the imdividual parameters on the *
C     screen.                                                              *
C     For IBM compatibel PCs and the MICROSOFT FORTRAN 5.0 Compiler the    *
C     following results are achieved:                                      *
C                                                                          *
C     NP =  21 : Number of points to be computed                           *
C     M  =  2 :                                                            *
C     K  =  3 : Order                                                      *
C     N  =  3 : N+1 is the number of DE BOOR points                        *
C     DP : DOUBLE PRECISION array DP(0:N+k-1,1:M):                         *
C                                                                          *
C        DP(0,I),I=1,2 = [0,0]                                             *
C        DP(1,I),I=1,2 = [2,3]                                             *
C        DP(2,I),I=1,2 = [5,1]                                             *
C        DP(3,I),I=1,2 = [3,0]                                             *
C                                                                          *
C                                                                          *
C--------------------------------------------------------------------------*
C                                                                          *
C     Required subroutines: CURVPG, KNOTVG, DEBOOR                         *
C                                                                          *
C--------------------------------------------------------------------------*
C[                                                                        ]*
C[  WANTED TO FIND: CLOSED B SPLINE                                       ]*
C[  ===============                                                       ]*
C[                                                                        ]*
C[                                                                        ]*
C[  SOLUTION: POINTS OF THE CLOSED B SPLINE                               ]*
C[  =========                                                             ]*
C[                                                                        ]*
C[                                                                        ]*
C[       I |         XP(I,1)      |         XP(I,2)      |                ]*
C[     ====+======================+======================|                ]*
C[       1 |  .10000000000000D+01 |  .15000000000000D+01 |                ]*
C[       2 |  .14200000423193D+01 |  .20000000596046D+01 |                ]*
C[       3 |  .18799999964237D+01 |  .22999999791384D+01 |                ]*
C[       4 |  .23800000619888D+01 |  .24000000000000D+01 |                ]*
C[       5 |  .29200000214577D+01 |  .22999999701977D+01 |                ]*
C[       6 |  .35000000000000D+01 |  .20000000000000D+01 |                ]*
C[       7 |  .40000000596046D+01 |  .16199999588728D+01 |                ]*
C[       8 |  .42999999791384D+01 |  .12800000083447D+01 |                ]*
C[       9 |  .44000000000000D+01 |  .97999996662140D+00 |                ]*
C[      10 |  .42999999701977D+01 |  .71999999761581D+00 |                ]*
C[      11 |  .40000000000000D+01 |  .50000000000000D+00 |                ]*
C[      12 |  .35799999576807D+01 |  .31999997973442D+00 |                ]*
C[      13 |  .31200000035763D+01 |  .18000000536442D+00 |                ]*
C[      14 |  .26199999380112D+01 |  .79999990463257D-01 |                ]*
C[      15 |  .20799999785423D+01 |  .20000003576278D-01 |                ]*
C[      16 |  .15000000000000D+01 |  .10649622220484D-28 |                ]*
C[      17 |  .99999994039536D+00 |  .60000001788139D-01 |                ]*
C[      18 |  .70000002086163D+00 |  .24000000715256D+00 |                ]*
C[      19 |  .60000000000000D+00 |  .54000004291535D+00 |                ]*
C[      20 |  .70000002980232D+00 |  .96000002861023D+00 |                ]*
C[      21 |  .10000000000000D+01 |  .15000000000000D+01 |                ]*
C                                                                          *
C                                                                          *
C***************************************************************************
C
      PARAMETER(NP=21,M=2,K=3,N=3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KV(N+2*(K-1)),DP(0:N+K-1,M),D(K,M),E(K,M),X(M)
     F  ,XP(1:NP,1:M)
      DATA (DP(0,MM),MM=1,2) / 0,0 /
      DATA (DP(1,MM),MM=1,2) / 2,3 /
      DATA (DP(2,MM),MM=1,2) / 5,1 /
      DATA (DP(3,MM),MM=1,2) / 3,0 /
      CALL CURVPG(NP,N,M,K,KV,DP,X,XP,D,E,IERR)
      IF (IERR .EQ. 0) THEN
        WRITE(*,890)
        WRITE(*,910)
        DO 10 I=1,NP
           WRITE(*,900)i,(XP(I,NN),NN=1,2)
10      CONTINUE
      ELSE
        WRITE(*,880) IERR
      ENDIF
      STOP
880   FORMAT(1X,'ERROR IERR =',I3)
890   FORMAT(1X,'C[',T76,']*',
     +    /,1X,'C[',2X,'WANTED TO FIND: CLOSED B SPLINE',T76,']*',
     +    /,1X,'C[',2X,'===============',T76,']*',
     +    /,1X,'C[',T76,']*',
     +    /,1X,'C[',T76,']*',
     +    /,1X,'C[  SOLUTION: POINTS OF THE CLOSED B SPLINE',T76,']*',
     +    /,1X,'C[  =========',T76,']*',
     +    /,1X,'C[',T76,']*')
900   FORMAT(1X,'C[',5X,I3,' |',1X,2(D20.14,1X,'|',1X),T76,']*')
910   FORMAT(1X,'C[',T76,']*',/,1X,'C[',7X,'I |',
     +    9X,'XP(I,1)',6X,'|',9X,'XP(I,2)',6X,'|',T76,']*',
     +    /,1X,'C[',5X,4('='),'+',22('='),'+',22('='),'|',T76,']*')
      END

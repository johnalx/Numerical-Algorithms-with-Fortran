      PROGRAM TEST
C
C*********************************************************************
C                                                                    *
C  Test program for the SUBROUTINE RKTRB                             *
C                                                                    *
C  Test examples:                                                    *
C  We solve two initial value problems using embedding formulas.     *
C                                                                    *
C  Test problem 1:                                                   *
C               Y'=COS(X)*3*Y^(2/3)           Y(0)=8    I=[0,20]     *
C                                                                    *
C               Each embedding formula solves the IVP using the      *
C               step size control HULL and CIVPS.                    *
C                                                                    *
C  Test problem 2:                                                   *
C               Y1'=Y2                        Y1(0)=1   I=[0,20]     *
C               Y2'=-Y1/(SQRT(Y1^2+Y3^2))^3   Y2(0)=0                *
C               Y3'=Y4                        Y3(0)=0                *
C               Y4'=-Y3/(SQRT(Y1^2+Y3^2))^3   Y4(0)=1                *
C                                                                    *
C               This problem is approached with the embedding formula*
C               RK8(7)13M and swep size control HULL. We compute the *
C               solution at X=5, X=10, X=15 and X=20 .               *
C                                                                    *
C  Each integration is continued until we reach the interval endpoint*
C  BETA. If IFEHL = -2 , i.e., if the maximally allowed number of    *
C  integration steps has been reached, we start RKTRB anew until we  *
C  reach BETA.                                                       *
C                                                                    *
C                                                                    *
C  Results:                                                          *
C                                                                    *
C[                                                                  ]*
C[ TEST 1:                                                          ]*
C[                                                                  ]*
C[ ABSERR= .10000000E-04    RELERR= .00000000E+00                   ]*
C[                                                                  ]*
C[ STEP CONTROL: 0                                                  ]*
C[                                                                  ]*
C[ METHOD I IFEHL I X     I APPROXIMATE VALUE I EXACT SOLUTION      ]*
C[ -------+-------+-------+-------------------+------------------   ]*
C[  0     I   0   I 20.00 I .247171948735E+02 I .247170687870E+02   ]*
C[  1     I   0   I 20.00 I .247182405223E+02 I .247170687870E+02   ]*
C[  2     I   0   I 20.00 I .247178976611E+02 I .247170687870E+02   ]*
C[  3     I   0   I 20.00 I .247172451880E+02 I .247170687870E+02   ]*
C[  4     I   0   I 20.00 I .247167894646E+02 I .247170687870E+02   ]*
C[  5     I   0   I 20.00 I .247174703964E+02 I .247170687870E+02   ]*
C[  6     I   0   I 20.00 I .247187273345E+02 I .247170687870E+02   ]*
C[  7     I   0   I 20.00 I .247170502638E+02 I .247170687870E+02   ]*
C[  8     I   0   I 20.00 I .247180503706E+02 I .247170687870E+02   ]*
C[  9     I   0   I 20.00 I .247172563100E+02 I .247170687870E+02   ]*
C[ 10     I   0   I 20.00 I .247190860505E+02 I .247170687870E+02   ]*
C[ 11     I   0   I 20.00 I .247170967601E+02 I .247170687870E+02   ]*
C[ 12     I   0   I 20.00 I .247170766655E+02 I .247170687870E+02   ]*
C[ 13     I   0   I 20.00 I .247175762869E+02 I .247170687870E+02   ]*
C[ 14     I   0   I 20.00 I .247173916666E+02 I .247170687870E+02   ]*
C[ 15     I   0   I 20.00 I .247170239069E+02 I .247170687870E+02   ]*
C[ 16     I   0   I 20.00 I .247170748700E+02 I .247170687870E+02   ]*
C[ 17     I   0   I 20.00 I .247169845490E+02 I .247170687870E+02   ]*
C[ 18     I   0   I 20.00 I .247169723416E+02 I .247170687870E+02   ]*
C[ 19     I   0   I 20.00 I .247171161728E+02 I .247170687870E+02   ]*
C[ 20     I   0   I 20.00 I .247170651334E+02 I .247170687870E+02   ]*
C[ 21     I   0   I 20.00 I .247170838516E+02 I .247170687870E+02   ]*
C[ 22     I   0   I 20.00 I .247170647267E+02 I .247170687870E+02   ]*
C[                                                                  ]*
C[ STEP CONTROL: 1                                                  ]*
C[                                                                  ]*
C[ METHOD I IFEHL I X     I APPROXIMATE VALUE I EXACT SOLUTION      ]*
C[ -------+-------+-------+-------------------+------------------   ]*
C[  0     I  -2   I  5.88 I .415290438787E+01 I .415290404059E+01   ]*
C[  0     I  -2   I 12.06 I .347899920452E+01 I .347899858857E+01   ]*
C[  0     I  -2   I 18.23 I .287335478745E+01 I .287335397647E+01   ]*
C[  0     I   0   I 20.00 I .247170720659E+02 I .247170687870E+02   ]*
C[  1     I   0   I 20.00 I .247172207487E+02 I .247170687870E+02   ]*
C[  2     I   0   I 20.00 I .247174720126E+02 I .247170687870E+02   ]*
C[  3     I   0   I 20.00 I .247171280422E+02 I .247170687870E+02   ]*
C[  4     I   0   I 20.00 I .247170169797E+02 I .247170687870E+02   ]*
C[  5     I   0   I 20.00 I .247172344202E+02 I .247170687870E+02   ]*
C[  6     I   0   I 20.00 I .247177932005E+02 I .247170687870E+02   ]*
C[  7     I   0   I 20.00 I .247170845840E+02 I .247170687870E+02   ]*
C[  8     I   0   I 20.00 I .247176421924E+02 I .247170687870E+02   ]*
C[  9     I   0   I 20.00 I .247172295387E+02 I .247170687870E+02   ]*
C[ 10     I  -2   I  2.84 I .121867747011E+02 I .121867548789E+02   ]*
C[ 10     I  -2   I  6.38 I .916202425104E+01 I .916188783637E+01   ]*
C[ 10     I  -2   I  9.64 I .568751769066E+01 I .568740411099E+01   ]*
C[ 10     I  -2   I 14.77 I .221234121931E+02 I .221228835315E+02   ]*
C[ 10     I  -2   I 18.48 I .440233726089E+01 I .440204842079E+01   ]*
C[ 10     I   0   I 20.00 I .247180040354E+02 I .247170687870E+02   ]*
C[ 11     I   0   I 20.00 I .247170883247E+02 I .247170687870E+02   ]*
C[ 12     I   0   I 20.00 I .247170736482E+02 I .247170687870E+02   ]*
C[ 13     I   0   I 20.00 I .247173746949E+02 I .247170687870E+02   ]*
C[ 14     I   0   I 20.00 I .247176148609E+02 I .247170687870E+02   ]*
C[ 15     I   0   I 20.00 I .247170234024E+02 I .247170687870E+02   ]*
C[ 16     I   0   I 20.00 I .247170701143E+02 I .247170687870E+02   ]*
C[ 17     I   0   I 20.00 I .247170363701E+02 I .247170687870E+02   ]*
C[ 18     I   0   I 20.00 I .247169824703E+02 I .247170687870E+02   ]*
C[ 19     I   0   I 20.00 I .247170882781E+02 I .247170687870E+02   ]*
C[ 20     I   0   I 20.00 I .247170863397E+02 I .247170687870E+02   ]*
C[ 21     I   0   I 20.00 I .247170788603E+02 I .247170687870E+02   ]*
C[ 22     I   0   I 20.00 I .247170961691E+02 I .247170687870E+02   ]*
C[                                                                  ]*
C[                                                                  ]*
C[ TEST 2:                                                          ]*
C[                                                                  ]*
C[ ABSERR= .00000000E+00    RELERR= .10000000E-05                   ]*
C[                                                                  ]*
C[ METHOD: RK8(7)13M   STEP CONTROL: HULL                           ]*
C[                                                                  ]*
C[ IFEHL=  0     X=   5.000                                         ]*
C[                                                                  ]*
C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
C[          ---------+----------------------+---------------------- ]*
C[          Y(1)     I  .28366237513788E+00 I  .28366218546323E+00  ]*
C[          Y(2)     I  .95892435839942E+00 I  .95892427466314E+00  ]*
C[          Y(3)     I -.95892378324852E+00 I -.95892427466314E+00  ]*
C[          Y(4)     I  .28366279894164E+00 I  .28366218546323E+00  ]*
C[                                                                  ]*
C[ IFEHL=  0     X=  10.000                                         ]*
C[                                                                  ]*
C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
C[          ---------+----------------------+---------------------- ]*
C[          Y(1)     I -.83906906063885E+00 I -.83907152907645E+00  ]*
C[          Y(2)     I  .54402455898843E+00 I  .54402111088937E+00  ]*
C[          Y(3)     I -.54402379555587E+00 I -.54402111088937E+00  ]*
C[          Y(4)     I -.83906961573268E+00 I -.83907152907645E+00  ]*
C[                                                                  ]*
C[ IFEHL=  0     X=  15.000                                         ]*
C[                                                                  ]*
C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
C[          ---------+----------------------+---------------------- ]*
C[          Y(1)     I -.75969379909921E+00 I -.75968791285882E+00  ]*
C[          Y(2)     I -.65028027791446E+00 I -.65028784015712E+00  ]*
C[          Y(3)     I  .65027965384846E+00 I  .65028784015712E+00  ]*
C[          Y(4)     I -.75969482356437E+00 I -.75968791285882E+00  ]*
C[                                                                  ]*
C[ IFEHL=  0     X=  20.000                                         ]*
C[                                                                  ]*
C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
C[          ---------+----------------------+---------------------- ]*
C[          Y(1)     I  .40806393616013E+00 I  .40808206181339E+00  ]*
C[          Y(2)     I -.91295369725738E+00 I -.91294525072763E+00  ]*
C[          Y(3)     I  .91295189103945E+00 I  .91294525072763E+00  ]*
C[          Y(4)     I  .40806470560004E+00 I  .40808206181339E+00  ]*
C                                                                    *
C*********************************************************************
C
C declarations
C
      EXTERNAL DGL1,DGL2
      DOUBLE PRECISION Y(4),WORK2(16,16),WORK1(4),YEX(4)
      INTEGER IWAHL(3),IWORK(2)
      DOUBLE PRECISION X,BETA,ABSERR,RELERR
C
C Test problem 1
C
      ABSERR=1.0D-05
      RELERR=0.0D+00
      N=1
      WRITE(*,1000) ABSERR,RELERR
C
C Loop choosing step size control
C
      DO 10 I=0,1
         IWAHL(2)=I
         WRITE(*,1100) IWAHL(2)
C
C loop choosing embedding formula
C
         DO 20 J=0,22
            IWAHL(1)=J
            IWAHL(3)=0
            X=0.0D+00
            BETA=20.D+00
            Y(1)=8.D+00
100         CALL RKTRB(X,BETA,N,DGL1,Y,ABSERR,RELERR,IWAHL,WORK1,
     +                 WORK2,IWORK,IFEHL)
            CALL EXAKT1(X,YEX(1))
            WRITE(*,1200) IWAHL(1),IFEHL,X,Y(1),YEX(1)
C
C if maximal number of integrations has been reached, call RKTRB
C anew
C
            IF(IFEHL.EQ.-2) THEN
              IWAHL(3)=1
              GOTO 100
            ENDIF
20       CONTINUE
10    CONTINUE
C
C Test problem 2
C
      ABSERR=0.0D+00
      RELERR=1.0D-06
      WRITE(*,2000) ABSERR,RELERR
      N=4
      IWAHL(1)=19
      IWAHL(2)=0
      IWAHL(3)=0
      X=0.0D+00
      Y(1)=1.0D+00
      Y(2)=0.0D+00
      Y(3)=0.0D+00
      Y(4)=1.0D+00
C
C loop to compute solution at  X=5, 10, 15 und 20
C
      DO 30 BETA=5.D+00,20.D+00+2.5D+00,5.D+00
200      CALL RKTRB(X,BETA,N,DGL2,Y,ABSERR,RELERR,IWAHL,WORK1,
     +              WORK2,IWORK,IFEHL)
         IWAHL(3)=1
         CALL EXAKT2(X,N,YEX)
         WRITE(*,2100) IFEHL,X
         DO 40 I=1,N
            WRITE(*,2200) I,Y(I),YEX(I)
40       CONTINUE
C
C if maximal number of integrations has been reached, call RKTRB
C anew
C
         IF(IFEHL.EQ.-2) THEN
           GOTO 200
         ELSEIF(IFEHL.GT.0.OR.IFEHL.LT.-2) THEN
           STOP
         ENDIF
30    CONTINUE
      STOP
C
C Format statements
C
1000  FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'TEST 1:',T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'ABSERR=',E14.8,'    RELERR=',E14.8,T70,']*')
1100  FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'STEP CONTROL: ',I1,T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'METHOD I',' IFEHL I X',5X,
     +       'I APPROXIMATE VALUE I EXACT SOLUTION',T70,']*',/,
     +       1X,'C[',1X,2(7('-'),'+'),7('-'),'+',19('-'),'+',18('-'),
     +       T70,']*')
1200  FORMAT(1X,'C[',1X,I2,5X,'I ',I3,3X,'I',F6.2,' I',E18.12,' I',
     +       E18.12,T70,']*')
2000  FORMAT(1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'TEST 2:',T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'ABSERR=',E14.8,'    RELERR=',E14.8,T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'METHOD: RK8(7)13M   STEP CONTROL: HULL',
     +       T70,']*')
2100  FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'IFEHL=',I3,5X,'X=',F8.3,T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[',10X,'SOLUTION I APPROXIMATE',10X,'I EXACT',
     +       T70,']*',/,
     +       1X,'C[',10X,9('-'),'+',22('-'),'+',22('-'),T70,']*')
2200  FORMAT(1X,'C[',10X,'Y(',I1,')',5X,'I ',E20.14,' I ',E20.14,
     +       T70,']*')
      END
C
C Test problems
C
      SUBROUTINE DGL1(X,Y,N,F)
      DOUBLE PRECISION Y(N),F(N),X
      F(1)=3.D+00*(ABS(Y(1)))**(2.D+00/3.D+00)*COS(X)
      RETURN
      END
C
      SUBROUTINE EXAKT1(X,YEXAKT)
      DOUBLE PRECISION X,YEXAKT
      YEXAKT=(2.D+00+SIN(X))**3.D+00
      RETURN
      END
C
C
      SUBROUTINE DGL2(X,Y,N,F)
      DOUBLE PRECISION Y(N),F(N),X
      F(1)=Y(2)
      F(2)=-Y(1)/(SQRT(Y(1)*Y(1)+Y(3)*Y(3)))**3.D+00
      F(3)=Y(4)
      F(4)=-Y(3)/(SQRT(Y(1)*Y(1)+Y(3)*Y(3)))**3.D+00
      RETURN
      END
C
      SUBROUTINE EXAKT2(X,N,YEXAKT)
      DOUBLE PRECISION X,YEXAKT(N)
      YEXAKT(1)=COS(X)
      YEXAKT(2)=-SIN(X)
      YEXAKT(3)=SIN(X)
      YEXAKT(4)=COS(X)
      RETURN
      END

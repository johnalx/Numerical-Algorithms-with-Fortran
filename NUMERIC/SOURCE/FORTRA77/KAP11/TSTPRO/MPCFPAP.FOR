      PROGRAM TEST
C
C*********************************************************************
C              Test program for the subroutines                      *
C              CFSPPA, CFSPPE, NCYFSY, PSPTAB.                       *
C--------------------------------------------------------------------*
C   required subroutiens    :  CFSPPA,CFSPPE,NCYFSY,NCYFSP,          *
C                              NCYFSS,PSPTAB                         *
C--------------------------------------------------------------------*
C                                                                    *
C   We compute a parametric, periodic cubic fitting spline for the   *
C   nodes:                                                           *
C                                                                    *
C          X            Y            WX          WY                  *
C      -------------------------------------------------             *
C        32.0          1.5          2.0         0.60                 *
C        25.0          2.5          0.2         0.08                 *
C        16.0          3.2          0.15        0.01                 *
C         5.0          3.0          0.15        0.01                 *
C        -6.0          3.1          0.40        0.10                 *
C       -18.0          1.8          0.10        0.20                 *
C       -16.0         -0.7          0.50        0.10                 *
C       -13.0         -2.2          0.50        0.10                 *
C        -6.0         -3.4          0.50        0.05                 *
C         4.0         -3.1          0.50        0.10                 *
C        21.0         -2.5          0.50        0.10                 *
C        33.0         -1.1          0.50        0.10                 *
C        32.0          1.5          2.00        0.60                 *
C                                                                    *
C      where WX(I) = weight for X(I), WY(I) = weight for Y(I).       *
C                                                                    *
C                                                                    *
C   The subroutine CFSPPA computes the parameters T(I), I=0, .., N.  *
C   These, the spline coefficients and the table of values for the   *
C   points XW, YW on the curve are put out.                          *  
C                                                                    *
C  Results:                                                          *
C                                                                    *
C[                                                                  ]*
C[ GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N=12),                ]*
C[ ======   AND WEIGHTS WX(I) AND WY(I)                             ]*
C[                                                                  ]*
C[  I        X(I)         Y(I)          WX(I)         WY(I)         ]*
C[  ---------------------------------------------------------       ]*
C[  0      32.0000       1.5000        2.0000         .6000         ]*
C[  1      25.0000       2.5000         .2000         .0800         ]*
C[  2      16.0000       3.2000         .1500         .0100         ]*
C[  3       5.0000       3.0000         .1500         .0100         ]*
C[  4      -6.0000       3.1000         .4000         .1000         ]*
C[  5     -18.0000       1.8000         .1000         .2000         ]*
C[  6     -16.0000       -.7000         .5000         .1000         ]*
C[  7     -13.0000      -2.2000         .5000         .1000         ]*
C[  8      -6.0000      -3.4000         .5000         .0500         ]*
C[  9       4.0000      -3.1000         .5000         .1000         ]*
C[ 10      21.0000      -2.5000         .5000         .1000         ]*
C[ 11      33.0000      -1.1000         .5000         .1000         ]*
C[ 12      32.0000       1.5000        2.0000         .6000         ]*
C[                                                                  ]*
C[                                                                  ]*
C[ TO FIND:   A) TABLE OF PARAMETER VALUES T(I),I=0(1)N             ]*
C[ ========                                                         ]*
C[            B) COEFFICIENTS OF A PARAMETRIC, PERIODIC, CUBIC      ]*
C[               APPROXIMATING SPLINE WITH COMPONENT FUNCTIONS      ]*
C[               SX AND SY FOR THE NODES X(I),Y(I).                 ]*
C[                                                                  ]*
C[            C) TABLE OF VALUES FOR SX(T),SY(T)                    ]*
C[                                                                  ]*
C[                                                                  ]*
C[ SOLUTION:                                                        ]*
C[ =========                                                        ]*
C[                                                                  ]*
C[ A) PARAMETER VALUES T:                                           ]*
C[    -------------------                                           ]*
C[                                                                  ]*
C[  I         T(I)                                                  ]*
C[ --------------------------                                       ]*
C[  0    .00000000000000E+00                                        ]*
C[  1    .70710678118655E+01                                        ]*
C[  2    .16098248988733E+02                                        ]*
C[  3    .27100067020314E+02                                        ]*
C[  4    .38100521556377E+02                                        ]*
C[  5    .50170732822134E+02                                        ]*
C[  6    .53372294940850E+02                                        ]*
C[  7    .56726396907100E+02                                        ]*
C[  8    .63828509268926E+02                                        ]*
C[  9    .73833008256881E+02                                        ]*
C[ 10    .90843593196853E+02                                        ]*
C[ 11    .10292498384528E+03                                        ]*
C[ 12    .10571066150072E+03                                        ]*
C[                                                                  ]*
C[                                                                  ]*
C[ B) SPLINE COEFFICIENTS:                                          ]*
C[    --------------------                                          ]*
C[                                                                  ]*
C[  I      AX(I)          BX(I)          CX(I)          DX(I)       ]*
C[ --------------------------------------------------------------   ]*
C[  0  .32003054E+02 -.59765569E+00 -.85889690E-01  .47521795E-02   ]*
C[  1  .25162655E+02 -.10994924E+01  .14919260E-01 -.66964087E-03   ]*
C[  2  .15960503E+02 -.99384179E+00 -.32156483E-02  .31777589E-03   ]*
C[  3  .50603842E+01 -.94920696E+00  .72726892E-02 -.11918288E-02   ]*
C[  4 -.60877766E+01 -.12218708E+01 -.32059288E-01  .46599443E-02   ]*
C[  5 -.17312179E+02  .40926295E-01  .13668025E+00 -.68037432E-02   ]*
C[  6 -.16003450E+02  .70689181E+00  .71332430E-01 -.65162578E-02   ]*
C[  7 -.13075855E+02  .96548060E+00  .57638508E-02 -.19501277E-03   ]*
C[  8 -.59980341E+01  .10178423E+01  .16088430E-02 -.35883522E-03   ]*
C[  9  .39866775E+01  .94228618E+00 -.91610569E-02  .75137555E-03   ]*
C[ 10  .21063081E+02  .12828704E+01  .29182956E-01 -.45053822E-02   ]*
C[ 11  .32876694E+02  .15195014E-01 -.13411089E+00  .57701343E-02   ]*
C[                                                                  ]*
C[                                                                  ]*
C[  I      AY(I)          BY(I)          CY(I)          DY(I)       ]*
C[ --------------------------------------------------------------   ]*
C[  0  .13412570E+01  .49561607E+00 -.63351536E-01  .29792086E-02   ]*
C[  1  .27315244E+01  .46571362E-01 -.15297651E-03 -.10778312E-03   ]*
C[  2  .30601785E+01  .17459729E-01 -.30719098E-02  .12525273E-03   ]*
C[  3  .30472373E+01 -.46516866E-02  .10621133E-02  .46523829E-04   ]*
C[  4  .31865239E+01  .35605318E-01  .25974631E-02 -.13955411E-02   ]*
C[  5  .15406403E+01 -.51163999E+00 -.47935966E-01  .72497834E-02   ]*
C[  6 -.35084160E+00 -.59564909E+00  .21695930E-01  .14304768E-02   ]*
C[  7 -.20506531E+01 -.40182977E+00  .36089825E-01 -.10586389E-02   ]*
C[  8 -.34633593E+01 -.49395034E-01  .13534107E-01 -.53064484E-03   ]*
C[  9 -.31342645E+01  .62072155E-01 -.23924005E-02  .40429902E-04   ]*
C[ 10 -.25716425E+01  .15776280E-01 -.32919168E-03  .12344709E-02   ]*
C[ 11 -.25222615E+00  .54837222E+00  .44413184E-01 -.12895093E-01   ]*
C[                                                                  ]*
C[ COEFFICIENTS FOUND.                                              ]*
C[                                                                  ]*
C[                                                                  ]*
C[ TABLE OF VALUES FOR:                                             ]*
C[   STARTING VALUE = T(0) =    .00                                 ]*
C[   FINAL VALUE    = T(N) = 105.71                                 ]*
C[                                                                  ]*
C[                                                                  ]*
C[ C) TABLE OF VALUES:                                              ]*
C[    ----------------                                              ]*
C[                                                                  ]*
C[   I          SX(T(I))                SY(T(I))                    ]*
C[ ---------------------------------------------------              ]*
C[   0     .32003053864485E+02     .13412569805887E+01              ]*
C[   1     .29026419423422E+02     .24332938819401E+01              ]*
C[   2     .25162654610217E+02     .27315243824472E+01              ]*
C[   3     .21971056534326E+02     .28673387149045E+01              ]*
C[   4     .18940162000859E+02     .29827632663872E+01              ]*
C[   5     .15960503329632E+02     .30601784905800E+01              ]*
C[   6     .12288240533373E+02     .30890718607117E+01              ]*
C[   7     .86235216008847E+01     .30724030283630E+01              ]*
C[   8     .50603841895773E+01     .30472373387416E+01              ]*
C[   9     .16188401531404E+01     .30467548901163E+01              ]*
C[  10    -.19796933793151E+01     .30885962694012E+01              ]*
C[  11    -.60877765978057E+01     .31865238973266E+01              ]*
C[  12    -.11219322846788E+02     .32809340756646E+01              ]*
C[  13    -.15567794261204E+02     .29140892897730E+01              ]*
C[  14    -.17312178746410E+02     .15406402638726E+01              ]*
C[  15    -.16003449825005E+02    -.35084160071688E+00              ]*
C[  16    -.13075854940430E+02    -.20506530579106E+01              ]*
C[  17    -.95834293703173E+01    -.30698850648783E+01              ]*
C[  18    -.59980341305443E+01    -.34633592926897E+01              ]*
C[  19    -.25991160599207E+01    -.34972492272745E+01              ]*
C[  20     .75573741907701E+00    -.33481908563177E+01              ]*
C[  21     .39866774706674E+01    -.31342644844405E+01              ]*
C[  22     .71159988321115E+01    -.29491863247286E+01              ]*
C[  23     .10210776297442E+02    -.28099372041646E+01              ]*
C[  24     .13448533444332E+02    -.27069649610614E+01              ]*
C[  25     .17006793850457E+02    -.26307174337319E+01              ]*
C[  26     .21063081093491E+02    -.25716424604890E+01              ]*
C[  27     .26408398466167E+02    -.24328234276140E+01              ]*
C[  28     .30934772334962E+02    -.18209339699404E+01              ]*
C[  29     .32876693801556E+02    -.25222614759190E+00              ]*
C[  30     .32003053864485E+02     .13412569805887E+01              ]*
C[                                                                  ]*
C[ TABLE FOUND.                                                     ]*
C                                                                    *
C*********************************************************************
C                                                                    *
C   Author      :  G�nter Palm                                       *
C   Date        :  5.22.1988                                         *
C   Source code :  FORTRAN 77                                        *
C                                                                    *
C*********************************************************************
C
      PARAMETER (N=12,MT=25)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION X(0:N),Y(0:N),T(0:N),WX(0:N),WY(0:N)
      DIMENSION AX(0:N),BX(0:N),CX(0:N),DX(0:N)
      DIMENSION AY(0:N),BY(0:N),CY(0:N),DY(0:N)
      DIMENSION XW(0:N+MT+2),YW(0:N+MT+2)
      DIMENSION ALPHA(2),BETA(2)
      DIMENSION HILF(14*N-10)
C
      DATA X / 32.0D0, 25.0 D0,16.0 D0, 5.0 D0,-6.0D0,-18.0D0,
     +        -16.0D0,-13.0 D0,-6.0 D0, 4.0 D0,21.0D0, 33.0D0,32.0D0/
      DATA Y /  1.5D0,  2.5 D0, 3.2 D0, 3.0 D0, 3.1D0,  1.8D0,
     +         -0.7D0, -2.2 D0,-3.4 D0,-3.1 D0,-2.5D0, -1.1D0, 1.5D0/
      DATA WX/  2.0D0,  0.2 D0, 0.15D0, 0.15D0, 0.4D0,  0.1D0,
     +          0.5D0,  0.5 D0, 0.5 D0, 0.5 D0, 0.5D0,  0.5D0, 2.0D0/
      DATA WY/  0.6D0,  0.08D0, 0.01D0, 0.01D0, 0.1D0,  0.2D0,
     +          0.1D0,  0.1 D0, 0.05D0, 0.1 D0, 0.1D0,  0.1D0, 0.6D0/
      DATA IW /2/ , JT /1/
      DATA IAB /4/
C
      CALL CFSPPA (N,X,Y,WX,WY,T,JT,IAB,ALPHA,BETA,IW,
     +             AX,BX,CX,DX,AY,BY,CY,DY,HILF,IFEHL)
      IF (IFEHL .NE. 0) THEN
         WRITE (*,1000) IFEHL
         STOP
      ENDIF
      WRITE (*,2000) N
      DO 100 I = 0,N,1
         WRITE (*,2100) I,X(I),Y(I),WX(I),WY(I)
  100 CONTINUE
      WRITE (*,2200)
      DO 150 I = 0,N,1
         WRITE (*,2300) I,T(I)
  150 CONTINUE
      WRITE (*,2400)
      DO 200 I = 0,N-1,1
         WRITE (*,2500) I,AX(I),BX(I),CX(I),DX(I)
  200 CONTINUE
      WRITE (*,2600)
      DO 250 I = 0,N-1,1
         WRITE (*,2500) I,AY(I),BY(I),CY(I),DY(I)
  250 CONTINUE
C
      WRITE (*,1050) T(0),T(N)
      TANF = T(0)
      TEND = T(N)
      CALL PSPTAB (N,MT,TANF,TEND,T,AX,BX,CX,DX,AY,BY,CY,DY,
     +             NT,XW,YW,IFEHL)
      IF (IFEHL .NE. 0) THEN
         WRITE (*,1300)
         STOP
      ENDIF
      WRITE (*,2700)
      DO 300 I = 0,NT,1
         WRITE (*,2800) I,XW(I),YW(I)
  300 CONTINUE
      WRITE (*,2900)
      STOP
C
C Format
C
 1000 FORMAT(1X,'C[ ERROR WITH COMPUTING COEFFICIENTS.',T70,']*',/,
     +       1X,'C[ ERROR PARAMETER =',I4,T70,']*')
 1050 FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[ COEFFICIENTS FOUND.',T70,']*',/,
     +       1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[ TABLE OF VALUES FOR:',T70,']*',/,
     +       1X,'C[   STARTING VALUE = T(0) = ',F6.2,T70,']*',/,
     +       1X,'C[   FINAL VALUE    = T(N) = ',F6.2,T70,']*',/,
     +       1X,'C[',T70,']*')
 1300 FORMAT(1X,'C[ ERRONEOUS INPUT.',T70,']*',/,
     +       1X,'C[ TANF EXCEEDS TEND.',T70,']*')
 2000 FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[ GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N=',I2,
     +          '),',T70,']*',/,
     +       1X,'C[ ======   AND WEIGHTS WX(I)',
     +          ' AND WY(I)',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[  I        X(I)         Y(I)          WX(I)',
     +          '         WY(I)',T70,']*',/,
     +       1X,'C[',2X,57('-'),T70,']*')
 2100 FORMAT(1X,'C[',1X,I2,4X,F9.4,4X,F9.4,4X,F10.4,4X,F10.4,T70,']*')
 2200 FORMAT(1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[ TO FIND:   A) TABLE OF PARAMETER VALUES ',
     +          'T(I),I=0(1)N',T70,']*',/,1X,'C[ ========',T70,']*',/,
     +       1X,'C[',12X,'B) COEFFICIENTS OF A PARAMETRIC, ',
     +          'PERIODIC, CUBIC',T70,']*',/,
     +       1X,'C[',15X,'APPROXIMATING SPLINE WITH COMPONENT ',
     +          'FUNCTIONS',T70,']*',/,
     +       1X,'C[',15X,'SX AND SY FOR THE NODES X(I),Y(I).',
     +          T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[',12X,'C) TABLE OF VALUES FOR SX(T),',
     +          'SY(T)',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[ SOLUTION:',T70,']*',/,
     +       1X,'C[ =========',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[ A) PARAMETER VALUES T:',T70,']*',/,
     +       1X,'C[    -------------------',T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[  I         T(I)',T70,']*',/,
     +       1X,'C[', 1X,26('-'),T70,']*')
 2300 FORMAT(1X,'C[',1X,I2,3X,E20.14,T70,']*')
 2400 FORMAT(1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[ B) SPLINE COEFFICIENTS:',T70,']*',/,
     +       1X,'C[    --------------------',T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[  I',6X,'AX(I)',10X,'BX(I)',10X,'CX(I)',10X,'DX(I)',
     +          T70,']*',/,
     +       1X,'C[',1X,62('-'),T70,']*')
 2500 FORMAT(1X,'C[',1X,I2,4(1X,E14.8),T70,']*')
 2600 FORMAT(1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[  I',6X,'AY(I)',10X,'BY(I)',10X,'CY(I)',10X,'DY(I)',
     +          T70,']*',/,
     +       1X,'C[',1X,62('-'),T70,']*')
 2700 FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[ C) TABLE OF VALUES:',T70,']*',/,
     +       1X,'C[    ----------------',T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[   I          SX(T(I))                SY(T(I))',
     +          T70,']*',/,1X,'C[',1X,51('-'),T70,']*')
 2800 FORMAT(1X,'C[',1X,I3,2(4X,E20.14),T70,']*')
 2900 FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[ TABLE FOUND.',T70,']*')
      END

      PROGRAM  TEST
C
C******************************************************************
C                                                                 *
C  Test program for the subroutine   FFT  with function values of *
C  the sum of the base functions                                  *
C                                                                 *
C         2 * EXP(I*(-3)*X) + I * EXP(I*7*X)                      *
C                                                                 *
C     ( I : imaginary unit (I**2 = -1), Period : 2*PI )           *
C  for 32 ( = 2**5) equidistant nodes in the interval [0, 2*PI).  *
C  The resulting discrete Fourier coefficients should all be zero,*
C  except                                                         *
C                                                                 *
C      RE ( F^(-3) ) = 2 ,  IM( F^(7) ) = 1  .                    *
C                                                                 *
C  After the inverse transformation (finding the functional       *
C  values), we display the difference to the original funactional *
C  values which should be close to zero throughout.               *
C  An IBM PC with Microsoft FORTRAN Compiler 5.0 gets the         *
C  following results:                                             *
C                                                                 *
C[   NODES:                                                      ]*
C[   ======                                                      ]*
C[   *    X    *     RE(F(X))     *     IM(F(X))     *           ]*
C[   *************************************************           ]*
C[   *  .00000 *  .2000000000E+01 *  .1000000000E+01 *           ]*
C[   *  .19635 *  .6821539442E+00 * -.9160501440E+00 *           ]*
C[   *  .39270 *  .3826834324E+00 * -.2771638598E+01 *           ]*
C[   *  .58905 *  .4412889683E+00 * -.2517140794E+01 *           ]*
C[   *  .78540 * -.7071067812E+00 * -.7071067812E+00 *           ]*
C[   *  .98175 * -.2517140794E+01 *  .4412889683E+00 *           ]*
C[   * 1.17810 * -.2771638598E+01 *  .3826834324E+00 *           ]*
C[   * 1.37445 * -.9160501440E+00 *  .6821539442E+00 *           ]*
C[   * 1.57080 *  .1000000000E+01 *  .2000000000E+01 *           ]*
C[   * 1.76715 *  .1306230788E+01 *  .2643724505E+01 *           ]*
C[   * 1.96350 *  .9238795325E+00 *  .1148050297E+01 *           ]*
C[   * 2.15984 *  .1406000328E+01 * -.1221650256E+01 *           ]*
C[   * 2.35619 *  .2121320344E+01 * -.2121320344E+01 *           ]*
C[   * 2.55254 *  .1221650256E+01 * -.1406000328E+01 *           ]*
C[   * 2.74889 * -.1148050297E+01 * -.9238795325E+00 *           ]*
C[   * 2.94524 * -.2643724505E+01 * -.1306230788E+01 *           ]*
C[   * 3.14159 * -.2000000000E+01 * -.1000000000E+01 *           ]*
C[   * 3.33794 * -.6821539442E+00 *  .9160501440E+00 *           ]*
C[   * 3.53429 * -.3826834324E+00 *  .2771638598E+01 *           ]*
C[   * 3.73064 * -.4412889683E+00 *  .2517140794E+01 *           ]*
C[   * 3.92699 *  .7071067812E+00 *  .7071067812E+00 *           ]*
C[   * 4.12334 *  .2517140794E+01 * -.4412889683E+00 *           ]*
C[   * 4.31969 *  .2771638598E+01 * -.3826834324E+00 *           ]*
C[   * 4.51604 *  .9160501440E+00 * -.6821539442E+00 *           ]*
C[   * 4.71239 * -.1000000000E+01 * -.2000000000E+01 *           ]*
C[   * 4.90874 * -.1306230788E+01 * -.2643724505E+01 *           ]*
C[   * 5.10509 * -.9238795325E+00 * -.1148050297E+01 *           ]*
C[   * 5.30144 * -.1406000328E+01 *  .1221650256E+01 *           ]*
C[   * 5.49779 * -.2121320344E+01 *  .2121320344E+01 *           ]*
C[   * 5.69414 * -.1221650256E+01 *  .1406000328E+01 *           ]*
C[   * 5.89049 *  .1148050297E+01 *  .9238795325E+00 *           ]*
C[   * 6.08684 *  .2643724505E+01 *  .1306230788E+01 *           ]*
C[                                                               ]*
C[   FOURIER COEFFICIENTS:                                       ]*
C[   =====================                                       ]*
C[   *  K  *     RE(C(K))     *     IM(C(K))     *               ]*
C[   *********************************************               ]*
C[   * -16 * -.2758210327E-15 * -.1561251128E-16 *               ]*
C[   * -15 * -.4608057773E-15 *  .4315681689E-16 *               ]*
C[   * -14 * -.2906118744E-15 * -.2688542524E-15 *               ]*
C[   * -13 *  .3896098871E-16 * -.5011104671E-15 *               ]*
C[   * -12 *  .2660991985E-16 * -.2663117816E-15 *               ]*
C[   * -11 *  .6497861939E-15 * -.9619254527E-16 *               ]*
C[   * -10 * -.6621028644E-16 *  .6866858617E-16 *               ]*
C[   *  -9 *  .1073936133E-15 *  .0000000000E+00 *               ]*
C[   *  -8 *  .8673617380E-16 * -.1734723476E-15 *               ]*
C[   *  -7 *  .3804896173E-15 * -.6497979704E-16 *               ]*
C[   *  -6 *  .2841488926E-15 *  .2672247160E-15 *               ]*
C[   *  -5 *  .1905113143E-15 *  .3657405656E-15 *               ]*
C[   *  -4 *  .6400891689E-17 *  .7587080969E-16 *               ]*
C[   *  -3 *  .2000000000E+01 *  .1221035263E-14 *               ]*
C[   *  -2 * -.4432072475E-15 *  .1160150479E-15 *               ]*
C[   *  -1 * -.2259469721E-15 * -.4092193168E-15 *               ]*
C[   *   0 * -.3642919300E-16 *  .2619432449E-15 *               ]*
C[   *   1 * -.2395944116E-15 * -.4847193093E-16 *               ]*
C[   *   2 * -.3549466923E-15 * -.3452327258E-15 *               ]*
C[   *   3 * -.3841138626E-16 * -.2331802483E-15 *               ]*
C[   *   4 *  .2891097528E-15 * -.5287733799E-16 *               ]*
C[   *   5 * -.4277415890E-15 * -.2820720402E-15 *               ]*
C[   *   6 *  .1061723824E-15 * -.4731691781E-15 *               ]*
C[   *   7 *  .4829476933E-15 *  .1000000000E+01 *               ]*
C[   *   8 *  .7285838599E-16 *  .9367506770E-16 *               ]*
C[   *   9 *  .3199105716E-15 * -.3988713471E-16 *               ]*
C[   *  10 *  .2920207351E-15 *  .2913511109E-15 *               ]*
C[   *  11 * -.3126306904E-16 *  .2575278474E-15 *               ]*
C[   *  12 * -.1972204740E-15 *  .1322960074E-15 *               ]*
C[   *  13 *  .3330669074E-15 * -.1094366568E-16 *               ]*
C[   *  14 *  .5630045633E-16 * -.4458136328E-16 *               ]*
C[   *  15 *  .3091933004E-16 *  .5757527705E-15 *               ]*
C[                                                               ]*
C[   INVERSE FOURIER TRANSFORM:                                  ]*
C[   ==========================                                  ]*
C[   *    X    *     RE(F(X))     *     IM(F(X))     *           ]*
C[   *************************************************           ]*
C[   *  .00000 *  .2000000000E+01 *  .1000000000E+01 *           ]*
C[   *  .19635 *  .6821539442E+00 * -.9160501440E+00 *           ]*
C[   *  .39270 *  .3826834324E+00 * -.2771638598E+01 *           ]*
C[   *  .58905 *  .4412889683E+00 * -.2517140794E+01 *           ]*
C[   *  .78540 * -.7071067812E+00 * -.7071067812E+00 *           ]*
C[   *  .98175 * -.2517140794E+01 *  .4412889683E+00 *           ]*
C[   * 1.17810 * -.2771638598E+01 *  .3826834324E+00 *           ]*
C[   * 1.37445 * -.9160501440E+00 *  .6821539442E+00 *           ]*
C[   * 1.57080 *  .1000000000E+01 *  .2000000000E+01 *           ]*
C[   * 1.76715 *  .1306230788E+01 *  .2643724505E+01 *           ]*
C[   * 1.96350 *  .9238795325E+00 *  .1148050297E+01 *           ]*
C[   * 2.15984 *  .1406000328E+01 * -.1221650256E+01 *           ]*
C[   * 2.35619 *  .2121320344E+01 * -.2121320344E+01 *           ]*
C[   * 2.55254 *  .1221650256E+01 * -.1406000328E+01 *           ]*
C[   * 2.74889 * -.1148050297E+01 * -.9238795325E+00 *           ]*
C[   * 2.94524 * -.2643724505E+01 * -.1306230788E+01 *           ]*
C[   * 3.14159 * -.2000000000E+01 * -.1000000000E+01 *           ]*
C[   * 3.33794 * -.6821539442E+00 *  .9160501440E+00 *           ]*
C[   * 3.53429 * -.3826834324E+00 *  .2771638598E+01 *           ]*
C[   * 3.73064 * -.4412889683E+00 *  .2517140794E+01 *           ]*
C[   * 3.92699 *  .7071067812E+00 *  .7071067812E+00 *           ]*
C[   * 4.12334 *  .2517140794E+01 * -.4412889683E+00 *           ]*
C[   * 4.31969 *  .2771638598E+01 * -.3826834324E+00 *           ]*
C[   * 4.51604 *  .9160501440E+00 * -.6821539442E+00 *           ]*
C[   * 4.71239 * -.1000000000E+01 * -.2000000000E+01 *           ]*
C[   * 4.90874 * -.1306230788E+01 * -.2643724505E+01 *           ]*
C[   * 5.10509 * -.9238795325E+00 * -.1148050297E+01 *           ]*
C[   * 5.30144 * -.1406000328E+01 *  .1221650256E+01 *           ]*
C[   * 5.49779 * -.2121320344E+01 *  .2121320344E+01 *           ]*
C[   * 5.69414 * -.1221650256E+01 *  .1406000328E+01 *           ]*
C[   * 5.89049 *  .1148050297E+01 *  .9238795325E+00 *           ]*
C[   * 6.08684 *  .2643724505E+01 *  .1306230788E+01 *           ]*
C[                                                               ]*
C[   ERROR TABLE:                                                ]*
C[   ============                                                ]*
C[   *   X    * ABS(RE(F(X)-F^(X))) * ABS(IM(F(X)-F^(X))) *      ]*
C[   ********************************************************    ]*
C[   *  .0000 *  .0000000000000E+00 *  .4440892098501E-15 *      ]*
C[   *  .1963 *  .4643095803669E-15 *  .7177418381854E-16 *      ]*
C[   *  .3927 *  .3017605696570E-15 *  .1987884683252E-15 *      ]*
C[   *  .5890 *  .2297966504583E-15 *  .1132882850030E-14 *      ]*
C[   *  .7854 *  .3252606517457E-17 *  .1100465205073E-16 *      ]*
C[   *  .9817 *  .7048398323328E-15 *  .5586351693732E-15 *      ]*
C[   * 1.1781 *  .5034492787936E-15 *  .1341700188451E-15 *      ]*
C[   * 1.3744 *  .6540314080247E-15 *  .3986611388229E-15 *      ]*
C[   * 1.5708 *  .7670730370335E-16 *  .1547698601223E-16 *      ]*
C[   * 1.7671 *  .7428140134241E-15 *  .8230720792424E-15 *      ]*
C[   * 1.9635 *  .4700016417725E-16 *  .5645982813218E-16 *      ]*
C[   * 2.1598 *  .1230569465771E-14 *  .6554002132675E-16 *      ]*
C[   * 2.3562 *  .5671461564272E-15 *  .6970877867996E-15 *      ]*
C[   * 2.5525 *  .6195673314668E-15 *  .1196525517555E-14 *      ]*
C[   * 2.7489 *  .7367966913668E-15 *  .2453549516335E-15 *      ]*
C[   * 2.9452 *  .1800968228716E-14 *  .1091669614948E-14 *      ]*
C[   * 3.1416 *  .3095397202446E-16 *  .1534146074067E-15 *      ]*
C[   * 3.3379 *  .4405113426809E-15 *  .2235624879665E-15 *      ]*
C[   * 3.5343 *  .3038476588391E-15 *  .4103705222858E-15 *      ]*
C[   * 3.7306 *  .2511012231476E-15 *  .8766858766718E-15 *      ]*
C[   * 3.9270 *  .1219727444046E-16 *  .2645453300865E-16 *      ]*
C[   * 4.1233 *  .1401006047286E-14 *  .5818913059730E-15 *      ]*
C[   * 4.3197 *  .6574059872866E-15 *  .5480641981914E-16 *      ]*
C[   * 4.5160 *  .7411063950025E-15 *  .3341511095600E-15 *      ]*
C[   * 4.7124 *  .8077306185017E-17 *  .4643095803669E-16 *      ]*
C[   * 4.9087 *  .5558975588876E-15 *  .1515280956266E-14 *      ]*
C[   * 5.1051 *  .1458794023079E-15 *  .1903316913798E-15 *      ]*
C[   * 5.3014 *  .1105235694632E-14 *  .6982261990807E-16 *      ]*
C[   * 5.4978 *  .4806810331714E-15 *  .7836613302725E-15 *      ]*
C[   * 5.6941 *  .9282939000821E-15 *  .1321859288694E-14 *      ]*
C[   * 5.8905 *  .5936006894358E-15 *  .2401507812055E-15 *      ]*
C[   * 6.0868 *  .1560004295881E-14 *  .8617509917458E-15 *      ]*
C                                                                 *
C******************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER  ITAU, M
      PARAMETER ( ITAU = 5 , M = 2**ITAU )
      DIMENSION  X(0:M-1), FUNKRE(0:M-1), FUNKIM(0:M-1)
C
C     Find real and imaginary parts of functional values
C
      PI = 4.0D0*ATAN(1.0D0)
      WRITE(*,900)
      DO 10 J = 0, M-1
        X(J) = DBLE(J)*2.0D0*PI/DBLE(M)
        FUNKRE(J) = 2.0D0*COS(3.0D0*X(J)) - SIN(7.0D0*X(J))
        FUNKIM(J) = -2.0D0*SIN(3.0D0*X(J)) + COS(7.0D0*X(J))
        WRITE(*,910) X(J), FUNKRE(J), FUNKIM(J)
  10  CONTINUE
C
C Compute Fourier coefficients
C
      CALL FFT (ITAU, FUNKRE, FUNKIM, 0)
      WRITE(*,920)
C
      DO 20 K = -M/2, -1
        WRITE (*, 940) K, FUNKRE(K+M), FUNKIM(K+M)
  20  CONTINUE
      DO 30 K = 0, M/2-1
        WRITE (*, 940) K, FUNKRE(K), FUNKIM(K)
  30  CONTINUE
C
C Inverse Fourier transformation
C
      CALL FFT (ITAU, FUNKRE, FUNKIM, 1)
C
      WRITE(*,930)
      DO 40 J = 0, M-1
        WRITE (*,910) X(J), FUNKRE(J), FUNKIM(J)
  40  CONTINUE
C
      WRITE(*,950)
      DO 50 J = 0, M-1
        WRITE (*, 960) X(J),
     #    DABS(2.0D0*COS(3.0D0*X(J))-SIN(7.0D0*X(J)) - FUNKRE(J)),
     #    DABS(-2.0D0*SIN(3.0D0*X(J))+COS(7.0D0*X(J)) - FUNKIM(J))
  50  CONTINUE
      STOP
C
C Format statements
C
  900 FORMAT(1X,'C[   NODES:',T67,']*',/,
     +       1X,'C[',3X,6('='),T67,']*',/,
     +       1X,'C[   *',4X,'X',4X,'*',5X,'RE(F(X))',5X,'*',
     +             5X,'IM(F(X))',5X,'*',T67,']*',/,
     +       1X,'C[',3X,49('*'),T67,']*')
  910 FORMAT(1X,'C[   * ',F7.5,' * ',2(E16.10,' * '),T67,']*')
  920 FORMAT(1X,'C[',T67,']*',/,
     +       1X,'C[   FOURIER COEFFICIENTS:',T67,']*',/,
     +       1X,'C[',3X,21('='),T67,']*',/,
     +       1X,'C[   *  K  *',5X,'RE(C(K))',5X,'*',5X,'IM(C(K))',
     +       5X,'*',T67,']*',/,1X,'C[',3X,45('*'),T67,']*')
  930 FORMAT(1X,'C[',T67,']*',/,
     +       1X,'C[   INVERSE FOURIER TRANSFORM:',T67,']*',/,
     +       1X,'C[',3X,26('='),T67,']*',/,
     +       1X,'C[   *',4X,'X',4X,'*',5X,'RE(F(X))',5X,'*',
     +             5X,'IM(F(X))',5X,'*',T67,']*',/,
     +       1X,'C[',3X,49('*'),T67,']*')
  940 FORMAT(1X,'C[   * ',I3,' * ',E16.10,' * ',E16.10,' *',T67,']*')
  950 FORMAT(1X,'C[',T67,']*',/,
     +       1X,'C[   ERROR TABLE:',T67,']*',/,
     +       1X,'C[',3X,12('='),T67,']*',/,
     +       1X,'C[   *',3X,'X',4X,'*',' ABS(RE(F(X)-F^(X))) *',
     +       ' ABS(IM(F(X)-F^(X))) *',T67,']*',/,
     +       1X,'C[',3X,56('*'),T67,']*')
  960 FORMAT(1X,'C[   * ',F6.4,' * ',2(E19.13,' * '),T67,']*')
      END

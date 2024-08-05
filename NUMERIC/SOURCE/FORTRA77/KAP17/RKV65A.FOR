C[BA*)
C[LE*)
      SUBROUTINE RKV65A (M,COEFF,QG)
C[IX{RKV65A}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKV6(5)9A (method of Verner).                          *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C M       : Dimension of the matrix COEFF depending on the chosen*
C           embedding formula                                    *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C COEFF   : 2-dimensional DOUBLE PRECISION array COEFF(1:16,1:M) *
C           with the coefficients for the embedding formula      *
C QG      : DOUBLE PRECISION value for the global error order    *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C H       : DOUBLE PRECISION auxiliary variable                  *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker KrÅger                                      *
C  Date     : 26.04.1993                                         *
C  Source   : FORTRAN 77                                         *
C  Sources  : J. H. Verner :                                     *
C             Some Runge-Kutta formular pairs                    *
C             page 496-511                                       *
C             SIAM J. NUMER. ANAL., Vol. 28, No. 2, 1991         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION COEFF(16,M),QG,H
C
C Initialize QG
C
      QG=5.0D0
C
C Initialize H
C
      H=DSQRT(10.0D0)
C
C Determine the matrix elements in COEFF
C
C         a - values (see scheme )
C
      COEFF(2,1)=1.0D0/8.0D0
      COEFF(3,1)=4.0D0/9.0D0-4.0D0/45.0D0*H
      COEFF(4,1)=2.0D0/3.0D0-2.0D0/15.0D0*H
      COEFF(5,1)=9.0D0/16.0D0
      COEFF(6,1)=1.0D0/2.0D0
      COEFF(7,1)=9.0D0/10.0D0
      COEFF(8,1)=1.0D0
      COEFF(9,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=1.0D0/8.0D0
      COEFF(3,2)=(-268.0D0+92.0D0*H)/405.0D0
      COEFF(4,2)=1.0D0/6.0D0-1.0D0/30.0D0*H
      COEFF(5,2)=11.547D3/32.768D3+405.0D0/16.384D3*H
      COEFF(6,2)=1.9662371D7/5.1149376D7+441.281D3/1.2787344D7*H
      COEFF(7,2)=-15.4207593D7/36.9412160D7-
     +           182.9424339D7/1154.4130000D7*H
      COEFF(8,2)=7059.4945601D7/2140.6013856D7+
     +           2147.3424323D7/2140.6013856D7*H
      COEFF(9,2)=31.0D0/324.0D0-37.0D0/4860.0D0*H
      COEFF(3,3)=(448.0D0-128.0D0*H)/405.0D0
      COEFF(4,4)=1.0D0/2.0D0-1.0D0/10.0D0*H
      COEFF(5,4)=-18.225D3/32.768D3-5103.0D0/16.384D3*H
      COEFF(6,4)=-3.786045D6/5.683264D6-2.52663D5/7.10408D5*H
      COEFF(7,4)=265.9895739D7/184.7060800D7+
     +           65.3855409D7/115.4413000D7*H
      COEFF(8,4)=-79.4525145D7/8.8090592D7-
     +           24.9156075D7/8.8090592D7*H
      COEFF(5,5)=12.555D3/16.384D3+2349.0D0/8192.0D0*H
      COEFF(6,5)=157.0556745D7/182.1486112D7+
     +           29.0041461D7/91.0743056D7*H
      COEFF(7,5)=-34949.2176711D7/59198.2986400D7-
     +           35978.4638379D7/147995.7466000D7*H
      COEFF(8,5)=86629.0968775D7/25409.7312624D7+
     +           25699.8959765D7/25409.731224D7*H
      COEFF(9,5)=(37.435D3-3235.0D0*H)/69.228D3
      COEFF(6,6)=(-4.1227072D7+1.374464D6*H)/51.2292969D7
      COEFF(7,6)=15392.0585664D7/9249.7341625D7+
     +           31106.6673408D7/46248.6708125D7*H
      COEFF(8,6)=(-1596419.6472448D7-503942.9245312D7*H)/
     +           128636.7645159D7
      COEFF(9,6)=-1.245184D6/1.090341D6+9.699328D6/1.6355115D7*H
      COEFF(7,7)=-1944.0D0/1625.0D0-6804.0D0/8125.0D0*H
      COEFF(8,7)=(17.017D3+5075.0D0*H)/1116.0D0
      COEFF(9,7)=71.0D0/54.0D0-74.0D0/135.0D0*H
      COEFF(8,8)=(42.875D3+16.625D3*H)/90.396D3
      COEFF(9,8)=625.0D0/486.0D0-250.0D0/729.0D0*H
      COEFF(9,9)=-23.0D0/21.0D0+37.0D0/105.0D0*H
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=31.0D0/324.0D0-37.0D0/4860.0D0*H
      COEFF(2,6)=(37.435D3-3235.0D0*H)/69.228D3
      COEFF(2,7)=-1.245184D6/1.090341D6+9.699328D6/1.6355115D7*H
      COEFF(2,8)=71.0D0/54.0D0-74.0D0/135.0D0*H
      COEFF(2,9)=625.0D0/486.0D0-250.0D0/729.0D0*H
      COEFF(3,4)=-23.0D0/21.0D0+37.0D0/105.0D0*H
C
C          A - values (see scheme)
C
      COEFF(1,1)=5.0D0/54.0D0-2.0D0/135.0D0*H
      COEFF(1,4)=(2390.0D0+2290.0D0*H)/17.307D3
      COEFF(1,5)=40.960D3/121.149D3+262.144D3/605.745D3*H
      COEFF(1,6)=2.0D0/27.0D0-64.0D0/135.0D0*H
      COEFF(1,8)=150.029D3/443.709D3-236.267D3/2.218545D6*H
      COEFF(1,9)=2411.0D0/126.774D3+1921.0D0/63.387D3*H
      RETURN
      END

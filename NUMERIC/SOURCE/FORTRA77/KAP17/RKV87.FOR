C[BA*)
C[LE*)
      SUBROUTINE RKV87 (M,COEFF,QG)
C[IX{RKV87}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKV8(7) (method of Verner).                            *
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
C             Explicit Runge-Kutta methods with estimates of     *
C             the local truncation error                         *
C             page 772-790                                       *
C             SIAM J. NUMER. ANAL., Vol. 15, No. 4, 1978         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION COEFF(16,M),QG
C
C Initialize QG
C
      QG=7.0D0
C
C Determine the matrix elements in COEFF
C
C         a - values (see scheme )
C
      COEFF(2,1)=1.0D0/4.0D0
      COEFF(3,1)=1.0D0/12.0D0
      COEFF(4,1)=1.0D0/8.0D0
      COEFF(5,1)=2.0D0/5.0D0
      COEFF(6,1)=1.0D0/2.0D0
      COEFF(7,1)=6.0D0/7.0D0
      COEFF(8,1)=1.0D0/7.0D0
      COEFF(9,1)=2.0D0/3.0D0
      COEFF(10,1)=2.0D0/7.0D0
      COEFF(11,1)=1.0D0
      COEFF(12,1)=1.0D0/3.0D0
      COEFF(13,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=1.0D0/4.0D0
      COEFF(3,2)=5.0D0/72.0D0
      COEFF(4,2)=1.0D0/32.0D0
      COEFF(5,2)=106.0D0/125.0D0
      COEFF(6,2)=1.0D0/48.0D0
      COEFF(7,2)=-1263.0D0/2401.0D0
      COEFF(8,2)=37.0D0/392.0D0
      COEFF(9,2)=1717.6D1/2551.5D1
      COEFF(10,2)=-2383.4D1/1800.75D2
      COEFF(11,2)=1273.3D1/7600.0D0
      COEFF(12,2)=-2706.1D1/2041.20D2
      COEFF(13,2)=1120.3D1/8680.0D0
      COEFF(3,3)=1.0D0/72.0D0
      COEFF(4,4)=3.0D0/32.0D0
      COEFF(5,4)=-408.0D0/125.0D0
      COEFF(5,5)=352.0D0/125.0D0
      COEFF(6,5)=8.0D0/33.0D0
      COEFF(7,5)=3993.6D1/2641.1D1
      COEFF(9,5)=-4710.4D1/2551.5D1
      COEFF(10,5)=-7782.4D1/1980.825D3
      COEFF(11,5)=-2003.2D1/5225.0D0
      COEFF(12,5)=4044.8D1/2806.65D2
      COEFF(13,5)=-3814.4D1/1193.5D1
      COEFF(6,6)=125.0D0/528.0D0
      COEFF(7,6)=-6412.5D1/2641.1D1
      COEFF(8,6)=1625.0D0/9408.0D0
      COEFF(9,6)=1325.0D0/504.0D0
      COEFF(10,6)=-6366.35D2/6338.64D2
      COEFF(11,6)=4564.85D2/8025.6D1
      COEFF(12,6)=-1353.775D3/1197.504D3
      COEFF(13,6)=2354.425D3/4583.04D2
      COEFF(7,7)=5520.0D0/2401.0D0
      COEFF(8,7)=-2.0D0/15.0D0
      COEFF(9,7)=-4179.2D1/2551.5D1
      COEFF(10,7)=2540.48D2/3001.25D2
      COEFF(11,7)=-4259.9D1/7125.0D0
      COEFF(12,7)=1766.2D1/2551.5D1
      COEFF(13,7)=-8404.6D1/1627.5D1
      COEFF(8,8)=61.0D0/6720.0D0
      COEFF(9,8)=2023.7D1/1458.00D2
      COEFF(10,8)=-183.0D0/7000.0D0
      COEFF(11,8)=3392.27D2/9120.00D2
      COEFF(12,8)=-7168.7D1/1166.400D3
      COEFF(13,8)=6733.09D2/1636.800D3
      COEFF(9,9)=4312.0D0/6075.0D0
      COEFF(10,9)=8.0D0/11.0D0
      COEFF(11,9)=-1029.0D0/4180.0D0
      COEFF(12,9)=98.0D0/225.0D0
      COEFF(13,9)=4704.0D0/8525.0D0
      COEFF(10,10)=-324.0D0/3773.0D0
      COEFF(11,10)=1701.0D0/1408.0D0
      COEFF(12,10)=1.0D0/16.0D0
      COEFF(13,10)=9477.0D0/1091.2D1
      COEFF(11,11)=5145.0D0/2432.0D0
      COEFF(12,11)=3773.0D0/1166.4D1
      COEFF(13,11)=-1029.0D0/992.0D0
      COEFF(13,13)=729.0D0/341.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=31.0D0/720.0D0
      COEFF(2,8)=16.0D0/75.0D0
      COEFF(2,9)=1680.7D1/7920.0D1
      COEFF(2,10)=1680.7D1/7920.0D1
      COEFF(2,11)=243.0D0/1760.0D0
      COEFF(3,4)=243.0D0/1760.0D0
      COEFF(3,5)=31.0D0/720.0D0
C
C          A - values (see scheme)
C
      COEFF(1,1)=13.0D0/288.0D0
      COEFF(1,6)=32.0D0/125.0D0
      COEFF(1,7)=3121.3D1/1440.00D2
      COEFF(1,8)=2401.0D0/1237.5D1
      COEFF(1,9)=1701.0D0/1408.0D1
      COEFF(1,10)=2401.0D0/1920.0D1
      COEFF(1,11)=19.0D0/450.0D0
      RETURN
      END

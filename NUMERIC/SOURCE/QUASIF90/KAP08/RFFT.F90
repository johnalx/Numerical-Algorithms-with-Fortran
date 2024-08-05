![             {Complex Discrete Fourier Transformation}*)              
      SUBROUTINE RFFT (ITAU, Y, ID) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     For ID = 0 this program determines the discrete Fourier    *      
!     coefficients                                               *      
!          A(0), ... , A(M/2) and B(1), ... , B(M/2-1)           *      
!     for M=2**ITAU given real function values Y(0), ... , Y(M-1)*      
!     of the associated discrete Fourier series                  *      
!                                                                *      
!          A(0) + (SUM K=1 to M/2-1)                             *      
!                     A(K)*COS(K*OMEGA*X) + B(K)*SIN(K*OMEGA*X)  *      
!               +  A(M/2)*COS(M/2*OMEGA*X)                       *      
!                                                                *      
!     Here OMEGA = 2*PI/L for the period L. For ID = 1, it per-  *      
!     forms the inverse transformation. Both transformations are *      
!     executed via a fast Fourier-transformation for the radix 2.*      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  ITAU: the number of function values is M=2**ITAU.             *      
!        ITAU has to be >= 2.                                    *      
!  Y   : real M-vector Y(0), Y(1), ... , Y(M-1);                 *      
!        depending on ID, Y has to be set up as follows:         *      
!         If ID=0 then Y contains the function values,           *      
!         If ID=1, then Y contains the discrete Fourier          *      
!                  coefficients labelled as follows:             *      
!                      Y(0)=A(0)                                 *      
!                      Y(K)=A((K+1)/2), K=1, 3, ... , M-1,       *      
!                      Y(K)=B(K/2)    , K=2, 4, ... , M-2,       *      
!                   in the following order:                      *      
!                      A(0), A(1), B(1), A(2), B(2), ...         *      
!  ID  : controls the direction of the performed transformation  *      
!           ID=0: compute the discrete Fourier coefficients      *      
!           ID=1: compute the function values                    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y   :   real M-vector Y(0), Y(1), ... , Y(M-1);               *      
!          depending on ID, the following values are returned:   *      
!           If ID=0: The vector contains the discrete Fourier    *      
!                    coefficients as follows:                    *      
!                      A(0)=Y(0)                                 *      
!                      A(K)=Y(2*K-1), K=1, ... , M/2,            *      
!                      B(K)=Y(2*K)  , K=1, ... , M/2-1,          *      
!                  in the following order:                       *      
!                      A(0), A(1), B(1), A(2), B(2), ...         *      
!           If ID=1: the function values are returned            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!                                                                *      
!  sources : Niederdrenk, K., see [NIED84].                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 02.14.1984                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER ITAU, SIGMA 
      DIMENSION Y (0:2**ITAU - 1) 
!                                                                       
      M = 2**ITAU 
      MD2 = M / 2 
      MD4 = MD2 / 2 
      FACTOR = 1.0D0 / DBLE (MD2) 
      PI = 4.0D0 * DATAN (1.0D0) 
      ARGMD2 = 2.0D0 * PI * FACTOR 
      ARGM = 0.5D0 * ARGMD2 
      VZ = DBLE (2 * ID-1) 
      IF (ID.EQ.1) FACTOR = 1.0D0 
!                                                                       
!*****************************************************************      
!     combining the real data for an FFT of halved length if     *      
!     ID = 1                                                     *      
!*****************************************************************      
!                                                                       
!****  URR, URI  : real and imaginary parts of      ****                
!****              the M-th roots of unity          ****                
!**** EPSR, EPSI : real and imaginary parts of      ****                
!****              (M-th root of unity)**K          ****                
                                                                        
      IF (ID.EQ.1) THEN 
         YSVE = Y (1) 
         Y (1) = Y (0) - Y (M - 1) 
         Y (0) = Y (0) + Y (M - 1) 
         URR = DCOS (ARGM) 
         URI = DSIN (ARGM) 
         EPSR = 1.0D0 
         EPSI = 0.0D0 
         DO 10 K = 1, MD4 - 1 
            STORE = EPSR 
            EPSR = STORE * URR - EPSI * URI 
            EPSI = STORE * URI + EPSI * URR 
            DMY1 = 0.5D0 * (EPSR * (YSVE-Y (M - 2 * K - 1) ) + EPSI *   &
            (Y (2 * K) + Y (M - 2 * K) ) )                              
            DMY2 = 0.5D0 * (EPSI * (YSVE-Y (M - 2 * K - 1) ) - EPSR *   &
            (Y (2 * K) + Y (M - 2 * K) ) )                              
            DMY3 = 0.5D0 * (YSVE+Y (M - 2 * K - 1) ) 
            DMY4 = 0.5D0 * (Y (2 * K) - Y (M - 2 * K) ) 
            YSVE = Y (2 * K + 1) 
            Y (2 * K) = DMY3 - DMY2 
            Y (2 * K + 1) = DMY1 - DMY4 
            Y (M - 2 * K) = DMY2 + DMY3 
            Y (M - 2 * K + 1) = DMY1 + DMY4 
   10    END DO 
         Y (MD2 + 1) = Y (MD2) 
         Y (MD2) = YSVE 
      ENDIF 
!                                                                       
!*****************************************************************      
!     relabelling with the bit-reversing function                *      
!     (while normalizing, if ID=0)                               *      
!*****************************************************************      
!                                                                       
      DO 30 J = 0, MD2 - 1 
         K = J 
         SIGMA = 0 
         DO 20 N = 1, ITAU - 1 
            KD2 = K / 2 
            SIGMA = 2 * SIGMA + K - 2 * KD2 
            K = KD2 
   20    END DO 
         IF (SIGMA.LT.J) GOTO 30 
         UR = Y (2 * J) 
         UI = Y (2 * J + 1) 
         Y (2 * J) = Y (2 * SIGMA) * FACTOR 
         Y (2 * J + 1) = Y (2 * SIGMA + 1) * FACTOR 
         Y (2 * SIGMA) = UR * FACTOR 
         Y (2 * SIGMA + 1) = UI * FACTOR 
   30 END DO 
!                                                                       
!*****************************************************************      
!      execution of the FFT of half length                       *      
!*****************************************************************      
!                                                                       
!****   MIN N   = 2**( ITAU-1 - N )                     ****            
!****  N MIN 1  = 2**( N - 1 )                          ****            
!****  N MIN 0  = 2**( N )                              ****            
!****   RR, RI  = real and imaginary parts of           ****            
!****             (M/2 -th root of unity)**(2**MINN)    ****            
!**** EPSR,EPSI = real and imaginary parts of           ****            
!****             (M/2 -th root of unity)**(L*2**MINN)  ****            
!                                                                       
      MINN = MD2 
      NMIN1 = 1 
      DO 130 N = 1, ITAU - 1 
         MINN = MINN / 2 
         NMIN0 = NMIN1 + NMIN1 
         ARG = ARGMD2 * DBLE (MINN) 
         RR = DCOS (ARG) 
         RI = VZ * DSIN (ARG) 
         EPSR = 1.0D0 
         EPSI = 0.0D0 
         DO 120 L = 0, NMIN1 - 1 
            DO 110 J = 0, MD2 - NMIN0, NMIN0 
               UR = Y (2 * (J + L) + NMIN0) * EPSR - Y (2 * (J + L)     &
               + NMIN0 + 1) * EPSI                                      
               UI = Y (2 * (J + L) + NMIN0) * EPSI + Y (2 * (J + L)     &
               + NMIN0 + 1) * EPSR                                      
               Y (2 * (J + L) + NMIN0) = Y (2 * (J + L) ) - UR 
               Y (2 * (J + L) + NMIN0 + 1) = Y (2 * (J + L) + 1)        &
               - UI                                                     
               Y (2 * (J + L) ) = Y (2 * (J + L) ) + UR 
               Y (2 * (J + L) + 1) = Y (2 * (J + L) + 1) + UI 
  110       END DO 
            STORE = EPSR 
            EPSR = STORE * RR - EPSI * RI 
            EPSI = STORE * RI + EPSI * RR 
  120    END DO 
         NMIN1 = NMIN0 
  130 END DO 
!                                                                       
!*****************************************************************      
!     separating the transformed data, if ID=0                   *      
!*****************************************************************      
!                                                                       
!****  URR, URI  : real and imaginary parts of    ****                  
!****              the M-th root of unity         ****                  
!**** EPSR, EPSI : real and imaginary parts of    ****                  
!****              (M-th root of unity)**K        ****                  
!                                                                       
      IF (ID.EQ.0) THEN 
         YSVE = Y (M - 1) 
         Y (M - 1) = 0.5D0 * (Y (0) - Y (1) ) 
         Y (0) = 0.5D0 * (Y (0) + Y (1) ) 
         URR = DCOS (ARGM) 
         URI = - DSIN (ARGM) 
         EPSR = 1.0D0 
         EPSI = 0.0D0 
         DO 150 K = 1, MD4 - 1 
            STORE = EPSR 
            EPSR = STORE * URR - EPSI * URI 
            EPSI = STORE * URI + EPSI * URR 
            DMY1 = 0.5D0 * (EPSI * (Y (2 * K) - Y (M - 2 * K) ) + EPSR *&
            (Y (2 * K + 1) + YSVE) )                                    
            DMY2 = 0.5D0 * (EPSR * (Y (2 * K) - Y (M - 2 * K) ) - EPSI *&
            (Y (2 * K + 1) + YSVE) )                                    
            DMY3 = 0.5D0 * (Y (2 * K) + Y (M - 2 * K) ) 
            DMY4 = 0.5D0 * (Y (2 * K + 1) - YSVE) 
            YSVE = Y (M - 2 * K - 1) 
            Y (2 * K - 1) = DMY1 + DMY3 
            Y (2 * K) = DMY2 - DMY4 
            Y (M - 2 * K - 1) = DMY3 - DMY1 
            Y (M - 2 * K) = DMY2 + DMY4 
  150    END DO 
         Y (MD2 - 1) = Y (MD2) 
         Y (MD2) = YSVE 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE RFFT                           

      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!  Test for SUBROUTINE CUTCHO:                                   *      
!  Read in a sparse positive definite matrix; reduce its band    *      
!  width via Cuthill-McKee; solve system via Cholesky method.    *      
!                                                                *      
!  Input file :  CUTCHO.MAT  Matrix element not zero             *      
!                            (structure see SUBROUTINE CUTCHO)   *      
!  Output files: CUTCHO.RS   some right hand sides used for      *      
!                            testing                             *      
!                CUTCHO.SOL  solutions from CUTCHO for thr right *      
!                            hand sides in CUTCHO.RS             *      
!                                                                *      
!  Compiling with Microsoft Fortran 5.0 and computing on a PC    *      
!  with coprocessor put the following output into CUTCHO.OUT.    *      
!                                                                *      
!----------------------------------------------------------------*      
!  Required subroutines:                                         *      
!                                                                *      
!  CUTHIL, CUTH1K, FNDROO, LVSTRU, SRTDEG, RDMTRX, BLDGPH,       *      
!  IBDWID, CUTPAK, PERMUT, CHOKYP, CHOKYS, MACHPD                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Elmar Pohl                                      *      
!  Date        : 11.17.1991                                      *      
!  Sourse code : Fortran 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      CHARACTER ( * ) FMAT, FSOL, FRS, FOUT 
      PARAMETER (FMAT = 'CUTCHO.MAT') 
      PARAMETER (FRS = 'CUTCHO.RS') 
      PARAMETER (FSOL = 'CUTCHO.SOL') 
      PARAMETER (FOUT = 'CUTCHO.OUT') 
      PARAMETER (IFRS = 7, IFSOL = 9, IFOUT = 10) 
!                                                                       
      PARAMETER (MAXELM = 50, MAXROW = 30, MAXAP = MAXROW * 5) 
      DIMENSION V (MAXELM), XEX (MAXROW), A (MAXROW, MAXROW) 
      INTEGER IC (MAXELM), IR (MAXROW) 
!                                                                       
      INTEGER NEIGHB (MAXELM), INB (MAXROW), LEVEL (MAXROW) 
      INTEGER ILV (MAXROW), IDEG (MAXROW), ICM (MAXROW) 
      INTEGER ICMREV (MAXROW) 
      LOGICAL MARK (MAXROW) 
      DIMENSION RSORG (MAXROW), RS (MAXROW), X (MAXROW), AP (MAXAP) 
      DIMENSION Z (MAXROW) 
!                                                                       
      NROW = 0 
      NV = 0 
      NLV = 0 
!                                                                       
!     Read matrix in order to display full matrix;                      
!     serves only as a demonstration. Usually one would not want to     
!     construct a full matrix (memory requirements).                    
!                                                                       
      CALL RDMTRX (FMAT, MAXELM, MAXROW, NROW, NV, V, IC, IR) 
!                                                                       
      OPEN (UNIT = IFOUT, FILE = FOUT) 
      WRITE (IFOUT, 1000) NROW, NV 
!                                                                       
!     construct full matrix and display                                 
!                                                                       
      DO 10 I = 1, NROW 
         DO 20 K = 1, NROW 
            A (I, K) = 0.D0 
   20    END DO 
   10 END DO 
!                                                                       
      DO 30 I = 1, NROW 
         DO 40 J = IR (I), IR (I + 1) - 1 
            A (I, IC (J) ) = V (J) 
   40    END DO 
   30 END DO 
      WRITE (IFOUT, 1100) 
!                                                                       
      DO 50 I = 1, NROW 
         WRITE (IFOUT, 1200) (A (I, K), K = 1, NROW) 
   50 END DO 
!                                                                       
!     Test: solve a number of equations with known solution.            
!     Create a file with suitable right hand sides.                     
!                                                                       
      NRS = 5 
      OPEN (UNIT = IFRS, FILE = FRS) 
      WRITE (IFRS, * ) NRS 
      DO 70 IRS = 1, NRS 
         DO 80 I = 1, NROW 
            XEX (I) = I + IRS - 1 
   80    END DO 
         DO 90 I = 1, NROW 
            RS (I) = 0.D0 
            DO 100 K = 1, NROW 
               RS (I) = RS (I) + A (I, K) * XEX (K) 
  100       END DO 
            WRITE (IFRS, * ) RS (I) 
   90    END DO 
   70 END DO 
      CLOSE (IFRS) 
!                                                                       
!     Solve linear system                                               
!                                                                       
      WRITE (IFOUT, 1300) 
      CALL CUTCHO (FMAT, FRS, FSOL, MAXELM, MAXROW, MAXAP, M, IFLAG, V, &
      IC, IR, NEIGHB, INB, LEVEL, ILV, IDEG, ICM, ICMREV, MARK, RSORG,  &
      RS, X, AP, Z)                                                     
      WRITE (IFOUT, 1400) M 
      IF (IFLAG.NE.1) THEN 
         WRITE (IFOUT, 1500) 
         STOP 
      ENDIF 
!                                                                       
!     Display solutions and compare with exact solutions                
!                                                                       
      OPEN (UNIT = IFSOL, FILE = FSOL) 
      OPEN (UNIT = IFRS, FILE = FRS) 
      READ (IFSOL, * ) NRS 
      READ (IFRS, * ) NRS 
      DO 110 IRS = 1, NRS 
         WRITE (IFOUT, 1600) IRS 
         WRITE (IFOUT, 1700) 
         WRITE (IFOUT, 1800) 
         DO 120 I = 1, NROW 
            XE = I + IRS - 1 
            READ (IFSOL, * ) XF 
            READ (IFRS, * ) RSI 
            WRITE (IFOUT, 1900) I, RSI, XF, XE, XE-XF 
  120    END DO 
  110 END DO 
      CLOSE (IFSOL) 
      CLOSE (IFRS) 
      CLOSE (IFOUT) 
!                                                                       
! Format statements                                                     
!                                                                       
 1000 FORMAT (1X,I3,' rows, ',I3,' elements') 
 1100 FORMAT (' Matrix:') 
 1200 FORMAT (1X,100G10.2) 
 1300 FORMAT (' CALL CUTCHO...') 
 1400 FORMAT (' Half band width after Cuthill-McKee: ',I4) 
 1500 FORMAT (' Matrix not positive definite') 
 1600 FORMAT (/,' Solution No',I3) 
 1700 FORMAT ('   i    R. H. S.        x(i)',                           &
     &        '            exact           error')                      
 1800 FORMAT (' ---------------------------',                           &
     &        '-----------------------------------')                    
 1900 FORMAT (1X,I3,3D16.8,D11.3) 
      END PROGRAM TEST                              

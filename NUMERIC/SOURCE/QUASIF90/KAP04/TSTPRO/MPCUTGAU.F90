      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!  Test for SUBROUTINE CUTGAU:                                   *      
!  read in a sparse nonsingular matrix, reduce its band width    *      
!  via  Cuthill-McKee, finally solve the system with Gauss       *      
!                                                                *      
!  Input file :  CUTGAU.MAT  Matrix elements not zero            *      
!                            (structure see SUBROUTINE RDMTRX)   *      
!  Output files: CUTGAU.RS   some right hand sides for testing   *      
!                CUTGAU.SOL  solutions from CUTGAU for the right *      
!                            hand sides in CUTGAU.RS             *      
!                                                                *      
!  When compiling with Microsoft Fortran 5.0 and computing on a  *      
!  PC with coprocessor the test program generated the output in  *      
!  CUTGAU.OUT.                                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!  Required subroutines:                                         *      
!                                                                *      
!  CUTHIL, CUTH1K, FNDROO, LVSTRU, SRTDEG, RDMTRX, BLDGPH,       *      
!  IBDWID, CUTPK2, PERMUT, BANDP,  BANDS,  MACHPD                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Elmar Pohl                                      *      
!  Date        : 11.17.1991                                      *      
!  Source code : Fortran 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      CHARACTER ( * ) FMAT, FSOL, FRS, FOUT 
      PARAMETER (FMAT = 'CUTGAU.MAT') 
      PARAMETER (FRS = 'CUTGAU.RS') 
      PARAMETER (FSOL = 'CUTGAU.SOL') 
      PARAMETER (FOUT = 'CUTGAU.OUT') 
      PARAMETER (IFRS = 7, IFSOL = 9, IFOUT = 10) 
!                                                                       
      PARAMETER (MAXELM = 100, MAXROW = 30, MAXAP = MAXROW * 16) 
      DIMENSION V (MAXELM), XEX (MAXROW), A (MAXROW, MAXROW) 
      INTEGER IC (MAXELM), IR (MAXROW) 
      INTEGER NEIGHB (MAXELM), INB (MAXROW), LEVEL (MAXROW) 
      INTEGER ILV (MAXROW), IDEG (MAXROW), ICM (MAXROW) 
      INTEGER ICMREV (MAXROW), IP (MAXROW) 
      LOGICAL MARK (MAXROW) 
      DIMENSION RSORG (MAXROW), RS (MAXROW), AP (MAXAP) 
!                                                                       
      NROW = 0 
      NV = 0 
      NLV = 0 
!                                                                       
!     read in matrix in order to display full matrix                    
!     (for demonstration purposes only. Usually avoid forming the full  
!     matrix (memory!) )                                                
!                                                                       
      CALL RDMTRX (FMAT, MAXELM, MAXROW, NROW, NV, V, IC, IR) 
!                                                                       
      OPEN (UNIT = IFOUT, FILE = FOUT) 
      WRITE (IFOUT, 1000) NROW, NV 
!                                                                       
!     form full matrix, display                                         
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
!     Test: solve a number of equations with given right hand sides     
!     Write out a file with suitable right hand sides                   
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
      CALL CUTGAU (FMAT, FRS, FSOL, MAXELM, MAXROW, MAXAP, M, IFLAG, V, &
      IC, IR, NEIGHB, INB, LEVEL, ILV, IDEG, ICM, ICMREV, MARK, RSORG,  &
      RS, AP, IP)                                                       
      WRITE (IFOUT, 1400) M 
      IF (IFLAG.NE.0) THEN 
         WRITE (IFOUT, 1500) 
         STOP 
      ENDIF 
!                                                                       
!     display and compare solution                                      
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
 1300 FORMAT (' CALL CUTGAU...') 
 1400 FORMAT (' Half band width after Cuthill-McKee: ',I4) 
 1500 FORMAT (' Matrix is numerically singular') 
 1600 FORMAT (/,' Solution No',I3) 
 1700 FORMAT ('   i    R. H. S.        x(i)',                           &
     &        '            exact           error')                      
 1800 FORMAT (' ---------------------------',                           &
     &        '-----------------------------------')                    
 1900 FORMAT (1X,I3,3D16.8,D11.3) 
      END PROGRAM TEST                              

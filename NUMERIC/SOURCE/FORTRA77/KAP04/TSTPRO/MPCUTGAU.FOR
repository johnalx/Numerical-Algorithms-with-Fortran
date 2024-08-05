      PROGRAM TEST
C
C*****************************************************************
C  Test for SUBROUTINE CUTGAU:                                   *
C  read in a sparse nonsingular matrix, reduce its band width    *
C  via  Cuthill-McKee, finally solve the system with Gauss       *
C                                                                *
C  Input file :  CUTGAU.MAT  Matrix elements not zero            *
C                            (structure see SUBROUTINE RDMTRX)   *
C  Output files: CUTGAU.RS   some right hand sides for testing   *
C                CUTGAU.SOL  solutions from CUTGAU for the right *
C                            hand sides in CUTGAU.RS             *
C                                                                *
C  When compiling with Microsoft Fortran 5.0 and computing on a  *
C  PC with coprocessor the test program generated the output in  *
C  CUTGAU.OUT.                                                   *
C                                                                *
C----------------------------------------------------------------*
C  Required subroutines:                                         *
C                                                                *
C  CUTHIL, CUTH1K, FNDROO, LVSTRU, SRTDEG, RDMTRX, BLDGPH,       *
C  IBDWID, CUTPK2, PERMUT, BANDP,  BANDS,  MACHPD                *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Elmar Pohl                                      *
C  Date        : 11.17.1991                                      *
C  Source code : Fortran 77                                      *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*(*) FMAT, FSOL, FRS, FOUT
      PARAMETER (FMAT='CUTGAU.MAT')
      PARAMETER (FRS ='CUTGAU.RS' )
      PARAMETER (FSOL='CUTGAU.SOL')
      PARAMETER (FOUT='CUTGAU.OUT')
      PARAMETER (IFRS=7, IFSOL=9, IFOUT=10)
C
      PARAMETER (MAXELM=100, MAXROW=30, MAXAP=MAXROW*16)
      DIMENSION V(MAXELM), XEX(MAXROW), A(MAXROW,MAXROW)
      INTEGER IC(MAXELM), IR(MAXROW)
      INTEGER NEIGHB(MAXELM), INB(MAXROW), LEVEL(MAXROW)
      INTEGER ILV(MAXROW), IDEG(MAXROW), ICM(MAXROW)
      INTEGER ICMREV(MAXROW), IP(MAXROW)
      LOGICAL MARK(MAXROW)
      DIMENSION RSORG(MAXROW), RS(MAXROW), AP(MAXAP)
C
      NROW = 0
      NV = 0
      NLV = 0
C
C     read in matrix in order to display full matrix
C     (for demonstration purposes only. Usually avoid forming the full
C     matrix (memory!) )
C
      CALL RDMTRX (FMAT, MAXELM, MAXROW, NROW, NV, V, IC, IR)
C
      OPEN (UNIT=IFOUT, FILE=FOUT)
      WRITE (IFOUT,1000) NROW,NV
C
C     form full matrix, display
C
      DO 10 I=1, NROW
         DO 20 K=1, NROW
            A(I,K) = 0.D0
   20    CONTINUE
   10 CONTINUE
C
      DO 30 I = 1, NROW
         DO 40 J = IR(I), IR(I + 1) - 1
            A(I, IC(J)) = V(J)
   40    CONTINUE
   30 CONTINUE
      WRITE (IFOUT,1100)
C
      DO 50 I = 1, NROW
         WRITE (IFOUT,1200) (A(I,K),K=1,NROW)
   50 CONTINUE
C
C     Test: solve a number of equations with given right hand sides
C     Write out a file with suitable right hand sides
C
      NRS = 5
      OPEN (UNIT=IFRS, FILE=FRS)
      WRITE (IFRS,*) NRS
      DO 70 IRS = 1, NRS
         DO 80 I = 1, NROW
               XEX(I) = I + IRS - 1
   80    CONTINUE
         DO 90 I = 1, NROW
            RS(I) = 0.D0
            DO 100 K = 1, NROW
               RS(I) = RS(I) + A(I, K) * XEX(K)
  100       CONTINUE
            WRITE (IFRS,*) RS(I)
   90    CONTINUE
   70 CONTINUE
      CLOSE (IFRS)
C
C     Solve linear system
C
      WRITE (IFOUT,1300)
      CALL CUTGAU (FMAT, FRS, FSOL, MAXELM, MAXROW, MAXAP, M,
     +             IFLAG, V, IC, IR, NEIGHB, INB, LEVEL, ILV,
     +             IDEG, ICM, ICMREV, MARK, RSORG, RS, AP, IP)
      WRITE (IFOUT,1400) M
      IF (IFLAG .NE. 0) THEN
         WRITE (IFOUT,1500)
         STOP
      END IF
C
C     display and compare solution
C
      OPEN (UNIT=IFSOL, FILE=FSOL)
      OPEN (UNIT=IFRS,   FILE=FRS  )
      READ (IFSOL,*) NRS
      READ (IFRS,*) NRS
      DO 110 IRS = 1, NRS
         WRITE (IFOUT,1600) IRS
         WRITE (IFOUT,1700)
         WRITE (IFOUT,1800)
         DO 120 I = 1, NROW
            XE = I + IRS - 1
            READ (IFSOL,*) XF
            READ (IFRS,*) RSI
            WRITE (IFOUT,1900) I, RSI, XF, XE, XE-XF
  120    CONTINUE
  110 CONTINUE
      CLOSE (IFSOL)
      CLOSE (IFRS)
      CLOSE (IFOUT)
C
C Format statements
C
1000  FORMAT (1X,I3,' rows, ',I3,' elements')
1100  FORMAT (' Matrix:')
1200  FORMAT (1X,100G10.2)
1300  FORMAT (' CALL CUTGAU...')
1400  FORMAT (' Half band width after Cuthill-McKee: ',I4)
1500  FORMAT (' Matrix is numerically singular')
1600  FORMAT (/,' Solution No',I3)
1700  FORMAT ('   i    R. H. S.        x(i)',
     +        '            exact           error')
1800  FORMAT (' ---------------------------',
     +        '-----------------------------------')
1900  FORMAT (1X,I3,3D16.8,D11.3)
      END

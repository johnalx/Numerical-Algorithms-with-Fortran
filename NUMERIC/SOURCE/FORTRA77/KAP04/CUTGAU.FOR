C[BA*)
C[LE*)
      SUBROUTINE CUTGAU (FMAT, FRS, FSOL, MAXELM, MAXROW, MAXAP,
C[IX{CUTGAU}*)
     +                   M, IFLAG, V, IC, IR, NEIGHB, INB, LEVEL,
     +                   ILV, IDEG, ICM, ICMREV, MARK, RSORG, RS,
     +                   AP, IP)
C
C*****************************************************************
C                                                                *
C  CUTGAU solves a linear system with a sparse and symmetric     *
C  system matrix using the Cuthill-McKee algorithm and Gau�      *
C  elimination with column pivot search.                         *
C[BE*)
C                                                                *
C  The nonzero elements of the system matrix are  read in        *
C  from the data file FMAT. The Cuthill-McKee method then trans- *
C  forms the matrix to one with band structure of minimal band   *
C  width.                                                        *
C  This matrix is then condensed and the resulting matrix is     *
C  factored using a standard Gau� algorithm for condensed        *
C  band matrices with column pivot search  (SUBROUTINE BANDP).   *
C  After BANDP, the right hand sides are read from the data set  *
C  FRS and solution vectors are found using BANDS. The solutions *
C  are stored in the data set FSOL.                              *
C                                                                *
C  For positive definite matrices the SUBROUTINE CUTCHO is       *
C  preferred over CUTGAU since CUTCHO uses the Cholesky method.  *
C  Hence the operations count is vey much reduced in CUTCHO and  *
C  the condensed matrix needs only about one third of the        *
C  storage that CUTGAU would use.                                *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FMAT   : CHAR*(*): Name of the input file of the nonzero      *
C           matrix elements                                      *
C           This data file is structured as follows:             *
C           - number of rows of the matrix                       *
C           - row wise:                                          *
C             - for every nonzero entry we store the tupel of    *
C               <column index> <entry>                           *
C               (indices start with 1)                           *
C             - a tupel with column index 0 designates the end of*
C               the row data                                     *
C  FRS    : CHAR*(*): Name of the input file for the right hand  *
C           sides in the form:                                   *
C           - number of right hand sides                         *
C           - for each right hand side:                          *
C             the entries appear in consecutive rows             *
C  FSOL   : CHAR*(*): Name of the output file for the solutions. *
C           It has the form:                                     *
C           - number of solutions (= number of right hand sides) *
C           - for each solution:                                 *
C             the entries appear in consecutive rows.            *
C  MAXELM : Dimensioning number for arrays, that will contain    *
C           matrix elements. MAXELM must be at least equal to the*
C           number of nonzero matrix elements.                   *
C  MAXROW : Dimensioning number for auxiliary arrays.            *
C           MAXROW must be at least as large as the number of    *
C           matrix rows.                                         *
C  MAXAP  : Upper index bound of vector AP                       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  M      : Number of upper or lower nonzero codiagonals of the  *
C           condensed matrix                                     *
C  IFLAG  : error parameter:                                     *
C               0: no error                                      *
C           .NE.0: Pivot element vanishes at the IFLAG-th        *
C                  elimination step; the matrix is numerically   *
C                  singular                                      *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  V      : vector V(1:MAXELM), Matrix in a linear list          *
C  IC     : vector IC(1:MAXELM), column indices of elements in V *
C  IR     : vector IR(1:MAXROW), starting indices of the rows of *
C           V                                                    *
C  NEIGHB : vector NEIGHB(1:MAXELM), incidence graph in a linear *
C           list                                                 *
C  INB    : vector INB(1:MAXROW), pointer for sublists in NEIGHB *
C  LEVEL  : vector LEVEL(1:MAXROW), level structure of a graph   *
C  ILV    : vector ILV(1:MAXROW), pointer for the levels in LEVEL*
C  IDEG   : vector IDEG(1:MAXROW), order of the nodes in the     *
C           graph                                                *
C  ICM    : vector ICM(1:MAXROW), Cuthill-McKee numbering        *
C  ICMREV : vector ICMREV(1:MAXROW), inverse permutation of ICM  *
C  MARK   : LOGICAL MARK(1:MAXROW), marks nodes                  *
C  RSORG  : vector RSORG(1:MAXROW), one right hand side          *
C  RS     : vector RS(1:MAXROW), RSORG after the Cuthill-McKee   *
C           permutation and performing BANDS                     *
C  AP     : vector AP(1:MAXROW*(3*MAXBB+1)), transformed matrix  *
C           in condensed form.                                   *
C           MAXBB is the number of upper codiagonals of the con- *
C           densed matrix. This number is found when performing  *
C           the Cuthill-McKee numbering scheme, hence MAXBB is   *
C           not known before the program is run. Estimates exist *
C           in the literature.                                   *
C  IP     : auxiliary vector IP(1:MAXROW) used in BANDP and BANDS*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines:                                         *
C                                                                *
C  RDMTRX  reads matrix from input file into arrays              *
C  BLDGPH  constructs the incidence graph of the matrix          *
C  CUTHIL  computes the Cuthill-McKee permutation                *
C  CUTH1K  Cuthill-McKee numbering for one component of the graph*
C  FNDROO  searches for starting nodes for the optimal level     *
C          structure                                             *
C  LVSTRU  constructs level structure of one component of the    *
C          graph                                                 *
C  SRTDEG  sorts nodes by their level                            *
C  IBDWID  computes half the band width of the condensed matrix  *
C  CUTPK2  forms condensed band matrix for SUBROUTINE BANDP      *
C  PERMUT  permutes the elements of a vector                     *
C  BANDP   decomposes a condensed band matrix using Gau�         *
C  BANDS   solves a linear system given in BANDP factorization   *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Elmar Pohl                                        *
C  Date      : 11.17.1991                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER IC(MAXELM), IR(MAXROW), NEIGHB(MAXELM), INB(MAXROW)
      INTEGER LEVEL(MAXROW), ILV(MAXROW), IP(MAXROW)
      INTEGER IDEG(MAXROW), ICM(MAXROW), ICMREV(MAXROW)
      LOGICAL MARK(MAXROW)
      DIMENSION V(MAXELM), RSORG(MAXROW), RS(MAXROW)
      DIMENSION AP(*)
      CHARACTER*(*) FMAT, FRS, FSOL
C
      PARAMETER (IFRS=8, IFSOL =9)
C
      NROW = 0
      NV = 0
      NLV = 0
C
C     Read matrix
C
      CALL RDMTRX (FMAT, MAXELM, MAXROW, NROW, NV, V, IC, IR)
C
C     Construct incidence graph
C
      CALL BLDGPH (NROW, IC, IR, NEIGHB, INB, IDEG)
C
C     Compute the Cuthill-McKee permutation
C
      CALL CUTHIL (NROW, NEIGHB, INB, IDEG, ICM,
     +             ICMREV, MARK, LEVEL, ILV)
C
C     Calculate half the band width of the condensed matrix
C
      M = IBDWID (NROW, NEIGHB, INB, ICM, ICMREV)
      if (nrow * (3 * m + 1) .gt. maxap) then
        write (*,*) 'CUTGAU: AP too small for condensed matrix!'
        stop
      end if
C
C     prevent the new permutation from increasing the bandwidth
C
      call ckbdwd(nrow, neighb, inb, m, icm, icmrev)
C
C     Condense the matrix for SUBROUTINE BANDP
C
      CALL CUTPK2 (NROW, M, V, IR, IC, ICMREV, AP)
C
C     Gau� factorization
C
      CALL BANDP (AP, NROW, 3*M+1, NROW, M, M, IP, ISIG, JFLAG)
      IFLAG = JFLAG
      IF (IFLAG .NE. 0) RETURN
C
C     Open data set FRS  in order to be able to read the right hand sides
C
      OPEN (UNIT=IFRS, FILE=FRS)
      READ (IFRS, *) NRS
C
C     Open data set FSOL for storing the solution vectors
C
      OPEN (UNIT=IFSOL , FILE=FSOL )
      WRITE (IFSOL , '(1X,I5)') NRS
C
C     Loop over all right hand sides
C
      DO 10 IRS = 1, NRS
C
C        Read the IRS-th right hand side
C
         READ (IFRS,*) (RSORG(I), I=1,NROW)
C
C        Permute right hand side
C
         CALL PERMUT (NROW, ICMREV, RSORG, RS)
C
C        Solve linear systems
C
         CALL BANDS (AP, NROW, 3*M+1, NROW, M, M, RS, IP)
C
C        Write solution into output data set
C
         DO 30 I = 1, NROW
            WRITE (IFSOL , '(1X,D17.10)') RS(ICMREV(I))
   30    CONTINUE
   10 CONTINUE
C
      CLOSE (IFRS)
      CLOSE (IFSOL )
C
      RETURN
      END

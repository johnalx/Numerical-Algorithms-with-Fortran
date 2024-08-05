      SUBROUTINE CUTGAU (FMAT, FRS, FSOL, MAXELM, MAXROW, MAXAP, M,     &
      IFLAG, V, IC, IR, NEIGHB, INB, LEVEL, ILV, IDEG, ICM, ICMREV,     &
      MARK, RSORG, RS, AP, IP)                                          
!                                                                       
!*****************************************************************      
!                                                                *      
!  CUTGAU solves a linear system with a sparse and symmetric     *      
!  system matrix using the Cuthill-McKee algorithm and Gauá      *      
!  elimination with column pivot search.                         *      
!                                                                *      
!  The nonzero elements of the system matrix are  read in        *      
!  from the data file FMAT. The Cuthill-McKee method then trans- *      
!  forms the matrix to one with band structure of minimal band   *      
!  width.                                                        *      
!  This matrix is then condensed and the resulting matrix is     *      
!  factored using a standard Gauá algorithm for condensed        *      
!  band matrices with column pivot search  (SUBROUTINE BANDP).   *      
!  After BANDP, the right hand sides are read from the data set  *      
!  FRS and solution vectors are found using BANDS. The solutions *      
!  are stored in the data set FSOL.                              *      
!                                                                *      
!  For positive definite matrices the SUBROUTINE CUTCHO is       *      
!  preferred over CUTGAU since CUTCHO uses the Cholesky method.  *      
!  Hence the operations count is vey much reduced in CUTCHO and  *      
!  the condensed matrix needs only about one third of the        *      
!  storage that CUTGAU would use.                                *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FMAT   : CHAR*(*): Name of the input file of the nonzero      *      
!           matrix elements                                      *      
!           This data file is structured as follows:             *      
!           - number of rows of the matrix                       *      
!           - row wise:                                          *      
!             - for every nonzero entry we store the tupel of    *      
!               <column index> <entry>                           *      
!               (indices start with 1)                           *      
!             - a tupel with column index 0 designates the end of*      
!               the row data                                     *      
!  FRS    : CHAR*(*): Name of the input file for the right hand  *      
!           sides in the form:                                   *      
!           - number of right hand sides                         *      
!           - for each right hand side:                          *      
!             the entries appear in consecutive rows             *      
!  FSOL   : CHAR*(*): Name of the output file for the solutions. *      
!           It has the form:                                     *      
!           - number of solutions (= number of right hand sides) *      
!           - for each solution:                                 *      
!             the entries appear in consecutive rows.            *      
!  MAXELM : Dimensioning number for arrays, that will contain    *      
!           matrix elements. MAXELM must be at least equal to the*      
!           number of nonzero matrix elements.                   *      
!  MAXROW : Dimensioning number for auxiliary arrays.            *      
!           MAXROW must be at least as large as the number of    *      
!           matrix rows.                                         *      
!  MAXAP  : Upper index bound of vector AP                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  M      : Number of upper or lower nonzero codiagonals of the  *      
!           condensed matrix                                     *      
!  IFLAG  : error parameter:                                     *      
!               0: no error                                      *      
!           .NE.0: Pivot element vanishes at the IFLAG-th        *      
!                  elimination step; the matrix is numerically   *      
!                  singular                                      *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  V      : vector V(1:MAXELM), Matrix in a linear list          *      
!  IC     : vector IC(1:MAXELM), column indices of elements in V *      
!  IR     : vector IR(1:MAXROW), starting indices of the rows of *      
!           V                                                    *      
!  NEIGHB : vector NEIGHB(1:MAXELM), incidence graph in a linear *      
!           list                                                 *      
!  INB    : vector INB(1:MAXROW), pointer for sublists in NEIGHB *      
!  LEVEL  : vector LEVEL(1:MAXROW), level structure of a graph   *      
!  ILV    : vector ILV(1:MAXROW), pointer for the levels in LEVEL*      
!  IDEG   : vector IDEG(1:MAXROW), order of the nodes in the     *      
!           graph                                                *      
!  ICM    : vector ICM(1:MAXROW), Cuthill-McKee numbering        *      
!  ICMREV : vector ICMREV(1:MAXROW), inverse permutation of ICM  *      
!  MARK   : LOGICAL MARK(1:MAXROW), marks nodes                  *      
!  RSORG  : vector RSORG(1:MAXROW), one right hand side          *      
!  RS     : vector RS(1:MAXROW), RSORG after the Cuthill-McKee   *      
!           permutation and performing BANDS                     *      
!  AP     : vector AP(1:MAXROW*(3*MAXBB+1)), transformed matrix  *      
!           in condensed form.                                   *      
!           MAXBB is the number of upper codiagonals of the con- *      
!           densed matrix. This number is found when performing  *      
!           the Cuthill-McKee numbering scheme, hence MAXBB is   *      
!           not known before the program is run. Estimates exist *      
!           in the literature.                                   *      
!  IP     : auxiliary vector IP(1:MAXROW) used in BANDP and BANDS*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:                                         *      
!                                                                *      
!  RDMTRX  reads matrix from input file into arrays              *      
!  BLDGPH  constructs the incidence graph of the matrix          *      
!  CUTHIL  computes the Cuthill-McKee permutation                *      
!  CUTH1K  Cuthill-McKee numbering for one component of the graph*      
!  FNDROO  searches for starting nodes for the optimal level     *      
!          structure                                             *      
!  LVSTRU  constructs level structure of one component of the    *      
!          graph                                                 *      
!  SRTDEG  sorts nodes by their level                            *      
!  IBDWID  computes half the band width of the condensed matrix  *      
!  CUTPK2  forms condensed band matrix for SUBROUTINE BANDP      *      
!  PERMUT  permutes the elements of a vector                     *      
!  BANDP   decomposes a condensed band matrix using Gauá         *      
!  BANDS   solves a linear system given in BANDP factorization   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      INTEGER IC (MAXELM), IR (MAXROW), NEIGHB (MAXELM), INB (MAXROW) 
      INTEGER LEVEL (MAXROW), ILV (MAXROW), IP (MAXROW) 
      INTEGER IDEG (MAXROW), ICM (MAXROW), ICMREV (MAXROW) 
      LOGICAL MARK (MAXROW) 
      DIMENSION V (MAXELM), RSORG (MAXROW), RS (MAXROW) 
      DIMENSION AP ( * ) 
      CHARACTER ( * ) FMAT, FRS, FSOL 
!                                                                       
      PARAMETER (IFRS = 8, IFSOL = 9) 
!                                                                       
      NROW = 0 
      NV = 0 
      NLV = 0 
!                                                                       
!     Read matrix                                                       
!                                                                       
      CALL RDMTRX (FMAT, MAXELM, MAXROW, NROW, NV, V, IC, IR) 
!                                                                       
!     Construct incidence graph                                         
!                                                                       
      CALL BLDGPH (NROW, IC, IR, NEIGHB, INB, IDEG) 
!                                                                       
!     Compute the Cuthill-McKee permutation                             
!                                                                       
      CALL CUTHIL (NROW, NEIGHB, INB, IDEG, ICM, ICMREV, MARK, LEVEL,   &
      ILV)                                                              
!                                                                       
!     Calculate half the band width of the condensed matrix             
!                                                                       
      M = IBDWID (NROW, NEIGHB, INB, ICM, ICMREV) 
      IF (nrow * (3 * m + 1) .gt.maxap) then 
         WRITE ( * , * ) 'CUTGAU: AP too small for condensed matrix!' 
         STOP 
      ENDIF 
!                                                                       
!     prevent the new permutation from increasing the bandwidth         
!                                                                       
      CALL ckbdwd (nrow, neighb, inb, m, icm, icmrev) 
!                                                                       
!     Condense the matrix for SUBROUTINE BANDP                          
!                                                                       
      CALL CUTPK2 (NROW, M, V, IR, IC, ICMREV, AP) 
!                                                                       
!     Gauá factorization                                                
!                                                                       
      CALL BANDP (AP, NROW, 3 * M + 1, NROW, M, M, IP, ISIG, JFLAG) 
      IFLAG = JFLAG 
      IF (IFLAG.NE.0) RETURN 
!                                                                       
!     Open data set FRS  in order to be able to read the right hand side
!                                                                       
      OPEN (UNIT = IFRS, FILE = FRS) 
      READ (IFRS, * ) NRS 
!                                                                       
!     Open data set FSOL for storing the solution vectors               
!                                                                       
      OPEN (UNIT = IFSOL, FILE = FSOL) 
      WRITE (IFSOL, '(1X,I5)') NRS 
!                                                                       
!     Loop over all right hand sides                                    
!                                                                       
      DO 10 IRS = 1, NRS 
!                                                                       
!        Read the IRS-th right hand side                                
!                                                                       
         READ (IFRS, * ) (RSORG (I), I = 1, NROW) 
!                                                                       
!        Permute right hand side                                        
!                                                                       
         CALL PERMUT (NROW, ICMREV, RSORG, RS) 
!                                                                       
!        Solve linear systems                                           
!                                                                       
         CALL BANDS (AP, NROW, 3 * M + 1, NROW, M, M, RS, IP) 
!                                                                       
!        Write solution into output data set                            
!                                                                       
         DO 30 I = 1, NROW 
            WRITE (IFSOL, '(1X,D17.10)') RS (ICMREV (I) ) 
   30    END DO 
   10 END DO 
!                                                                       
      CLOSE (IFRS) 
      CLOSE (IFSOL) 
!                                                                       
      RETURN 
      END SUBROUTINE CUTGAU                         

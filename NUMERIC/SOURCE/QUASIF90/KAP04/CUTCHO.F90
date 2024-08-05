![  {The Cuthill--McKee Algorithm}                                      
![  {The Algorithm of Cuthill--McKee for Sparse Symmetric               
![   Matrices}*)                                                        
      SUBROUTINE CUTCHO (FMAT, FRS, FSOL, MAXELM, MAXROW, MAXAP, M,     &
      IFLAG, V, IC, IR, NEIGHB, INB, LEVEL, ILV, IDEG, ICM, ICMREV,     &
      MARK, RSORG, RS, X, AP, Z)                                        
!                                                                       
!*****************************************************************      
!                                                                *      
!  CUTCHO solves a linear system with a sparse and symmetric     *      
!  positive definite system matrix using the Cuthill-McKee       *      
!  algorithm and the Cholesky decomposition.                     *      
!                                                                *      
!  The nonzero elements of the system matrix are  read in        *      
!  from the data file FMAT. The Cuthill-McKee method then trans- *      
!  forms the matrix to one with band structure of minimal band   *      
!  width.                                                        *      
!  The upper half band of this symmetric matrix is then condensed*      
!  and the resulting matrix is factored using a Cholesky method  *      
!  for condensed symmetric band matrices (SUBROUTINE CHOBDZ).    *      
!  After CHOBDZ, the right hand sides are read from the data set *      
!  FRS and solution vectors are found using CHOBDL. The solutions*      
!  are stored in the data set FSOL.                              *      
!                                                                *      
!  If the given system matrix is not positive definite, its      *      
!  Cholesky decomposition will be obtainable and thus the        *      
!  SUBROUTINE CUTCHO cannot be used. In this case we recommend   *      
!  using CUTGAU which used a version of the Gauá algorithm for   *      
!  condensed matrices with partial pivoting that, however, uses  *      
!  more operations and three times the storage of CUTCHO.        *      
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
!               1: no error                                      *      
!           .NE.1: the matrix is numerically singular or not     *      
!                  positive definite                             *      
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
!           permutation                                          *      
!  X      : vector X(1:MAXROW), the solution still to be permuted*      
!  AP     : vector AP(1:MAXROW*(M+1)), upper band of the         *      
!           transformed matrix in condensed form.                *      
!           M is the number of upper codiagonals of the condensed*      
!           matrix. This number is found when performing the     *      
!           Cuthill-McKee numbering scheme, hence M is not known *      
!           before the program is run. Estimates exist in the    *      
!           literature.                                          *      
!  Z      : auxiliary vector IP(1:MAXROW) used in CHOBDZ         *      
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
!  CUTPAK  forms condensed band matrix for SUBROUTINE CHOBDZ     *      
!  PERMUT  permutes the elements of a vector                     *      
!  CHOBDZ  decomposes a condensed band matrix using Cholesky     *      
!  CHOBDL  solves a linear system given in CHOBDZ factorization  *      
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
      INTEGER LEVEL (MAXROW), ILV (MAXROW) 
      INTEGER IDEG (MAXROW), ICM (MAXROW), ICMREV (MAXROW) 
      LOGICAL MARK (MAXROW) 
      DIMENSION V (MAXELM), RSORG (MAXROW), RS (MAXROW) 
      DIMENSION X (MAXROW), Z (MAXROW), AP ( * ) 
      CHARACTER ( * ) FMAT, FRS, FSOL 
!                                                                       
      PARAMETER (IFRS = 8, IFSOL = 9) 
!                                                                       
      NROW = 0 
      NV = 0 
      NLV = 0 
!                                                                       
!     Read in matrix                                                    
!                                                                       
      CALL RDMTRX (FMAT, MAXELM, MAXROW, NROW, NV, V, IC, IR) 
!                                                                       
!     Construct the incidence gragh of the matrix                       
!                                                                       
      CALL BLDGPH (NROW, IC, IR, NEIGHB, INB, IDEG) 
!                                                                       
!     Find the CUTHILL-MCKEE permutation                                
!                                                                       
      CALL CUTHIL (NROW, NEIGHB, INB, IDEG, ICM, ICMREV, MARK, LEVEL,   &
      ILV)                                                              
!                                                                       
!     Find half the bandwidth of the transformed matrix                 
!                                                                       
      M = IBDWID (NROW, NEIGHB, INB, ICM, ICMREV) 
      IF (nrow * m.gt.maxap) then 
         WRITE ( * , * ) 'CUTCHO: AP too small for condensed matrix!' 
         STOP 
      ENDIF 
!                                                                       
!     prevent the new permutation from increasing the band width        
!                                                                       
      CALL ckbdwd (nrow, neighb, inb, m, icm, icmrev) 
!                                                                       
!     Condense the matrix for Cholesky                                  
!                                                                       
      CALL CUTPAK (NROW, M, V, IR, IC, ICMREV, AP) 
!                                                                       
!     Cholesky decomposition                                            
!                                                                       
      CALL CHOBDZ (NROW, M, AP, JFLAG, Z) 
      IFLAG = JFLAG 
      IF (IFLAG.NE.1) RETURN 
!                                                                       
!     open the data file FRS, prepare to read in right hand sides       
!                                                                       
      OPEN (UNIT = IFRS, FILE = FRS) 
      READ (IFRS, * ) NRS 
!                                                                       
!     open the data set FSOL for the solution vectors                   
!                                                                       
      OPEN (UNIT = IFSOL, FILE = FSOL) 
      WRITE (IFSOL, '(1X,I5)') NRS 
!                                                                       
!     Loop for all right hand sides                                     
!                                                                       
      DO 10 IRS = 1, NRS 
!                                                                       
!        Read in the  IRS-th right hand side                            
!                                                                       
         READ (IFRS, * ) (RSORG (I), I = 1, NROW) 
!                                                                       
!        Permute right hand side according to Cuthill-McKee permutaion  
!                                                                       
         CALL PERMUT (NROW, ICMREV, RSORG, RS) 
!                                                                       
!        Solve linear system                                            
!                                                                       
         CALL CHOBDL (NROW, M, AP, RS, X, Z) 
!                                                                       
!        Write solution onto output data file                           
!                                                                       
         DO 30 I = 1, NROW 
            WRITE (IFSOL, '(1X,D17.10)') X (ICMREV (I) ) 
   30    END DO 
   10 END DO 
!                                                                       
      CLOSE (IFRS) 
      CLOSE (IFSOL) 
!                                                                       
      RETURN 
      END SUBROUTINE CUTCHO                         
!                                                                       
!                                                                       
      SUBROUTINE CUTHIL (NNODES, NEIGHB, INB, IDEG, ICM, ICMREV, MARK,  &
      LEVEL, ILV)                                                       
!                                                                       
!*****************************************************************      
!                                                                *      
!  CUTHIL computes the Cuthill-McKee numbering of a graph.       *      
!  This Cuthill-McKee numbering is used to solve linear syatems  *      
!  with sparse and symmetric system matrices and saves storage   *      
!  space and computation al time. When using the Cuthill-McKee   *      
!  permutation on the graph of a symmetric matrix, the matrix is *      
!  transformed into a symmetric band matrix with a generally     *      
!  reduced bandwidth.                                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NNODES : number of nodes of the graph. The graph is determined*      
!           by the following two vectors:                        *      
!  NEIGHB : vector NEIGHB(1:*), the list of adjacent nodes.      *      
!           For I=1, ..., NNODES, the vector NEIGHB contains the *      
!           numbers of adjacent nodes of node I in positions     *      
!           NEIGHB(K) where K=INB(I), ..., INB(I+1)-1.           *      
!  INB    : vector INB(1:NNODES+1) containing the indices for    *      
!           NEIGHB.                                              *      
!           INB(NNODES+1) must be equal to the number of entries *      
!           in NEIGHB plus 1.                                    *      
!  IDEG   : vector IDEG(1:NNODES), containing the degree of every*      
!           node, i.e., the number of its neighbors.             *      
!                                                                *      
!  The vectors  NEIGHB, INB and IDEG can be formed from A by     *      
!  using SUBROUTINE BLDGPH.                                      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ICM    : vector ICM(1:NNODES) containing the permutations of  *      
!           the nodes according to Cuthill-McKee. For I=1, ...,  *      
!           NNODES,  ICM(I)  describes the original index of the *      
!           node, while I is its Cuthill-McKee number            *      
!  ICMREV : vector ICMREV(1:NNODES), the inverse permutation of  *      
!           ICM: For I=1(1)NNODES, ICMREV(I) denotes the Cuthill-*      
!           McKee number of the node originally numbered I.      *      
!                                                                *      
!  REMARK:                                                       *      
!        One of the permutation vectors ICM or ICMREV is clearly *      
!        redundant. One can form the transformed matrix          *      
!        completely from ICM. However, if one wants to condense  *      
!        the transformed matrix, one would have to conduct       *      
!        expensive searches inside ICM, unless one has its       *      
!        inverse ICMREV available.                               *      
!                                                                *      
!                                                                *      
!  AUXILIARY VECTORS:                                            *      
!  ==================                                            *      
!  MARK   : LOGICAL MARK(1:NNODES) for labelling nodes           *      
!  LEVEL  : vector LEVEL(1:NNODES) for denoting the level        *      
!           structure of one component of the graph              *      
!  ILV    : vector ILV(1:NNODES) containing the starting indices *      
!           of levels in the vector LEVEL                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:                                         *      
!                                                                *      
!  FNDROO   searches for starting nodes for an optimal level     *      
!           structure                                            *      
!  CUTH1K   Cuthill-McKee numbering of one component of the graph*      
!  LVSTRU   constructs the level structure of graph component    *      
!  SRTDEG   sorts nodes according to degree                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NEIGHB ( * ), INB (1:NNODES + 1), IDEG (1:NNODES) 
      INTEGER ICM (1:NNODES), ICMREV (1:NNODES) 
      INTEGER LEVEL (1:NNODES), ILV (1:NNODES) 
      LOGICAL MARK (1:NNODES) 
                                                                        
      DO 10 I = 1, NNODES 
         MARK (I) = .FALSE. 
         ICM (I) = 0 
   10 END DO 
      NFOUND = 0 
      DO 20 I = 1, NNODES 
         IF (.NOT.MARK (I) ) THEN 
!                                                                       
!           Start of a new component in the graph                       
!                                                                       
            IROOT = I 
!                                                                       
!           Search a starting node that will give a level structure of  
!           maximal length                                              
!                                                                       
            CALL FNDROO (IROOT, NNODES, NEIGHB, INB, IDEG, MARK, NLV,   &
            LEVEL, ILV, LNODES)                                         
!                                                                       
!           Cuthill-McKee numbering of this component                   
!                                                                       
            CALL CUTH1K (IROOT, NFOUND+1, NNODES, NEIGHB, INB, IDEG,    &
            MARK, ICM)                                                  
            NFOUND = NFOUND+LNODES 
         ENDIF 
   20 END DO 
!                                                                       
!     All components of the graph have been numbered.                   
!     Form the inverse ICMREV of ICM                                    
!                                                                       
      DO 30 I = 1, NNODES 
         ICMREV (ICM (I) ) = I 
   30 END DO 
      RETURN 
      END SUBROUTINE CUTHIL                         
!                                                                       
!                                                                       
      SUBROUTINE CUTH1K (IROOT, ISTART, NNODES, NEIGHB, INB, IDEG, MARK,&
      ICM)                                                              
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for the Cuthill-McKee algorithm.                   *      
!  It determines the component induced by the node IROOT together*      
!  with its  Cuthill-McKee numbering.                            *      
!                                                                *      
!  In order to find the Cuthill-McKee numbering of a graph, one  *      
!  calls CUTHIL and not CUTH1K. CUTHIL in turn calls CUTH1K      *      
!  repeatedly until all components of the graph have been found. *      
!  Hence  CUTHIL will work equally well for disconnected graphs. *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IROOT  : Number of the starting node                          *      
!           (output of SUBROUTINE FNDROO())                      *      
!  ISTART : Starting index for the Cuthill-McKee numbering of the*      
!           new component                                        *      
!  NNODES : Number of nodes of the graph. The graph is stored in *      
!           the two vectors below:                               *      
!  NEIGHB : vector NEIGHB(1:*), the list of adjacent nodes.      *      
!           For I=1, ..., NNODES, the vector NEIGHB contains the *      
!           numbers of adjacent nodes of node I in positions     *      
!           NEIGHB(K) where K=INB(I), ..., INB(I+1)-1.           *      
!  INB    : vector INB(1:NNODES+1) containing the indices for    *      
!           NEIGHB.                                              *      
!           INB(NNODES+1) must be equal to the number of entries *      
!           in NEIGHB plus 1.                                    *      
!  IDEG   : vector IDEG(1:NNODES), containing the degree of every*      
!           node, i.e., the number of its neighbors.             *      
!  MARK   : LOGICAL MARK(1:NNODES), node markers.                *      
!           If MARK(I)=.FALSE., then the node labelled I can be  *      
!           used in the new component. Otherwise it belongs to   *      
!           another componentn from an earlier call of CUTH1K.   *      
!                                                                *      
!  The vectors  NEIGHB, INB and IDEG can be formed by using      *      
!  SUBROUTINE BLDGPH.                                            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ICM    : vector ICM(1:NNODES) containing the permutations of  *      
!           the nodes according to Cuthill-McKee. For I=1, ...,  *      
!           NNODES,  ICM(I)  describes the original index of the *      
!           node, while I is its Cuthill-McKee number            *      
!           The entries for I=ISTART, ..., ISTART+NNEW-1  are    *      
!           adjusted when calling CUTH1K, where NNEW denotes the *      
!           number of nodes of the new component.                *      
!  MARK   : same as on input, except that the indices of nodes in*      
!           the new component are now labelled .TRUE. .          *      
!                                                                *      
!  By successively calling  CUTH1K  for increasing values of     *      
!  ISTART, all positions in ICM and MARK are set.                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:                                         *      
!                                                                *      
!  SRTDEG   sorts nodes according to degree                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NEIGHB ( * ), INB (1:NNODES + 1), IDEG (1:NNODES) 
      INTEGER ICM (1:NNODES) 
      LOGICAL MARK (1:NNODES) 
!                                                                       
!     The method is similar to the algorithm for finding the            
!     level structure. This level structure is formed in the            
!     vector ICM. In addition we also make lists of nodes               
!     (ordered by degree) that are appended to existing level sets.     
!                                                                       
      ICM (ISTART) = IROOT 
      MARK (IROOT) = .TRUE. 
                                                                        
      NEWEND = ISTART 
!                                                                       
!     NEWEND denotes the final node of a component                      
!     while it is being formed                                          
!                                                                       
      LEVEND = ISTART - 1 
   10 CONTINUE 
      LEVBEG = LEVEND+1 
      LEVEND = NEWEND 
!                                                                       
!        LEVBEG denotes the start of the last level in ICM(),           
!        LEVEND marks the end.                                          
!                                                                       
!        Find nodes of the next level:                                  
!        search for nodes adjacent to lower level nodes                 
!        and enter them in ICM in case they have not been               
!        marked .TRUE. before.                                          
!                                                                       
      DO 20 I = LEVBEG, LEVEND 
!                                                                       
!           Find a list of unmarked neighbors of the                    
!           node originally labelled by ICM(I).                         
!                                                                       
         NEWBEG = NEWEND+1 
!                                                                       
!           ICM is the starting index of this list                      
!                                                                       
         DO 30 J = INB (ICM (I) ), INB (ICM (I) + 1) - 1 
            IF (.NOT.MARK (NEIGHB (J) ) ) THEN 
               NEWEND = NEWEND+1 
               ICM (NEWEND) = NEIGHB (J) 
               MARK (NEIGHB (J) ) = .TRUE. 
            ENDIF 
   30    END DO 
!                                                                       
!           Sort ICM(NEWBEG), ..., ICM(NEWEND) by increasing degree     
!                                                                       
         CALL SRTDEG (ICM, IDEG, NEWBEG, NEWEND) 
   20 END DO 
!                                                                       
!     Stay inside this loop as long as new nodes are being found.       
!                                                                       
      IF (NEWEND.GT.LEVEND) GOTO 10 
      RETURN 
      END SUBROUTINE CUTH1K                         
!                                                                       
!                                                                       
      SUBROUTINE BLDGPH (NROW, IC, IR, NEIGHB, INB, IDEG) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for the Cuthill-McKee method.                      *      
!  Form the graph of a symmetric matrix.                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NROW   : number of rows                                       *      
!  IC     : vector IC(1:NV) with the column indices of the       *      
!           nonzero matrix entries                               *      
!           (NV = number of nonzero matrix entries)              *      
!  IR     : vector IR(1:NROW+1) with the indices for the         *      
!           beginnings of rows. IR(NROW+1) must equal NV+1.      *      
!                                                                *      
!  All inputs are available as outputs of SUBROUTINE RDMTRX.     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NEIGHB : vector NEIGB(1:NV-NROW) with the indices of adjacent *      
!           nodes.                                               *      
!           Each row I corresponds to node I. If I.NE.K, then    *      
!           node I is adjacent to node K provided A(I,K).NE.0.   *      
!           Since the system matrix is assumed symmetric, if I   *      
!           is adjacent to K, so is K to I.                      *      
!           For I=1, ..., NROW, the vector NEIGHB contains the   *      
!           indices of adjacent nodes for node I from index      *      
!           INB(I) to index  INB(I+1)-1).                        *      
!  INB    : vector INB(1:NROW+1) with indices for the vector     *      
!           NEIGHB.                                              *      
!           For I=1, ..., NROW, INB(I) denotes the starting      *      
!           index of the list of adjacent nodes for I in NEIGHB. *      
!           This list extends to index INB(I+1)-1. We always have*      
!           INB(NROW+1) = NV-NROW+1.                             *      
!  IDEG   : vector IDEG(1:NROW) which specifies the degree of    *      
!           each node I=1, ..., NROW, i.e., the number of nodes  *      
!           adjacent to node I.                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER IC ( * ), IR (1:NROW + 1), NEIGHB ( * ) 
      INTEGER INB (1:NROW + 1), IDEG (1:NROW) 
!                                                                       
      N = 0 
      DO 10 I = 1, NROW 
         INB (I) = N + 1 
         DO 20 K = IR (I), IR (I + 1) - 1 
            IF (IC (K) .NE.I) THEN 
               N = N + 1 
               NEIGHB (N) = IC (K) 
            ENDIF 
   20    END DO 
         IDEG (I) = N + 1 - INB (I) 
   10 END DO 
      INB (NROW + 1) = N + 1 
      RETURN 
      END SUBROUTINE BLDGPH                         
!                                                                       
!                                                                       
      SUBROUTINE CUTPAK (N, M, V, IR, IC, ICMREV, AP) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for solving a linear system via  Cuthill-McKee     *      
!  and Cholesky.                                                 *      
!  CUTPAK performs the Cuthill-McKee permutation specified in ICM*      
!  for the matrix that is given in  V, IR and IC and condenses   *      
!  the upper band of the resulting symmetric matrix in the array *      
!  AP.                                                           *      
!  This array AP can subsequently be used to solve linear systems*      
!  using the  SUBROUTINE CHOBND/CHOBDZ/CHOBDL.                   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : number of rows of the matrix                         *      
!  M      : number of upper codiagonals after  Cuthill-McKee     *      
!           permutation                                          *      
!  V      : vector V(1:*), which contains the nonzero elements of*      
!           the matrix row after row                             *      
!  IC     : vector IC(1:*) of column indices for each entry in V *      
!  IR     : vector IR(1:N+1), containing the starting index for  *      
!           each row in V or IC. We must have                    *      
!               IR(N+1) =  1 + number of entries in V .          *      
!  ICMREV : vector ICMREV(1:N), the inverse permutation of the   *      
!           Cuthill-McKee numbering ICM. For I=1, ..., N,        *      
!           ICMREV(I) denotes the Cuthill-McKee number of the    *      
!           node with the original index I.                      *      
!                                                                *      
!  N,V,IR and IC are outputs of SUBROUTINE RDMTRX.               *      
!  M is an output of FUNCTION IBDWID.                            *      
!  ICMREV is an output of SUBROUTINE CUTHIL.                     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  AP     : array AP(1:N,1:M+1), the condensed matrix            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
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
      DIMENSION V ( * ), AP (1:N, * ) 
      INTEGER IC ( * ), IR ( * ), ICMREV ( * ) 
!                                                                       
      DO 10 I = 1, N 
         DO 20 K = 1, M + 1 
            AP (I, K) = 0.0D0 
   20    END DO 
   10 END DO 
!                                                                       
      DO 30 I = 1, N 
         IREV = ICMREV (I) 
         DO 40 K = IR (I), IR (I + 1) - 1 
            KREV = ICMREV (IC (K) ) 
            IF (KREV.GE.IREV) AP (IREV, KREV - IREV + 1) = V (K) 
   40    END DO 
   30 END DO 
      RETURN 
      END SUBROUTINE CUTPAK                         
!                                                                       
!                                                                       
      SUBROUTINE CUTPK2 (N, M, V, IR, IC, ICMREV, AP) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for solving a linear system via  Cuthill-McKee     *      
!  and Gauá.                                                     *      
!  CUTPK2 performs the Cuthill-McKee permutation specified in ICM*      
!  for the matrix that is given in  V, IR and IC and condenses   *      
!  the upper band of the resulting symmetric matrix in the array *      
!  AP.                                                           *      
!  This array AP can subsequently be used to solve linear systems*      
!  using the  SUBROUTINE BAND/BANDP/BANDS.                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : number of rows of the matrix                         *      
!  M      : number of upper codiagonals after  Cuthill-McKee     *      
!           permutation                                          *      
!  V      : vector V(1:*), which contains the nonzero elements of*      
!           the matrix row after row                             *      
!  IC     : vector IC(1:*) of column indices for each entry in V *      
!  IR     : vector IR(1:N+1), containing the starting index for  *      
!           each row in V or IC. We must have                    *      
!               IR(N+1) =  1 + number of entries in V .          *      
!  ICMREV : vector ICMREV(1:N), the inverse permutation of the   *      
!           Cuthill-McKee numbering ICM. For I=1, ..., N,        *      
!           ICMREV(I) denotes the Cuthill-McKee number of the    *      
!           node with the original index I.                      *      
!                                                                *      
!  N,V,IR and IC are outputs of SUBROUTINE RDMTRX.               *      
!  M is an output of FUNCTION IBDWID.                            *      
!  ICMREV is an output of SUBROUTINE CUTHIL.                     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  AP     : array AP(1:N,1:M+1), the condensed matrix            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
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
      DIMENSION V ( * ), AP (1:N, * ) 
      INTEGER IC ( * ), IR ( * ), ICMREV ( * ) 
!                                                                       
      DO 10 I = 1, N 
         DO 20 K = 1, 2 * M + 1 
            AP (I, K) = 0.0D0 
   20    END DO 
   10 END DO 
!                                                                       
      DO 30 I = 1, N 
         IREV = ICMREV (I) 
         DO 40 K = IR (I), IR (I + 1) - 1 
            AP (IREV, M + 1 + ICMREV (IC (K) ) - IREV) = V (K) 
   40    END DO 
   30 END DO 
      RETURN 
      END SUBROUTINE CUTPK2                         
!                                                                       
!                                                                       
      SUBROUTINE FNDROO (IROOT, NNODES, NEIGHB, INB, IDEG, MARK, NLV,   &
      LEVEL, ILV, LNODES)                                               
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for the Cuthill-McKee algorithm.                   *      
!  It constructs the level structure of the component of the     *      
!  graph of IROOT and attempts to choose a starting node so that *      
!  the resulting structure will have as many levels as possible. *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IROOT  : Number of the node that defines the component        *      
!  NNODES : Number of nodes of the graph. The graph is defined by*      
!           the following two vectors:                           *      
!  NEIGHB : vector NEIGHB(1:*) with the lists of adjacent nodes  *      
!           For I=1, ..., NNODES, the vector NEIGHB contains the *      
!           indices of the nodes that are adjacent to node I in  *      
!           NEIGHB(K) for  K=INB(I), ...,  INB(I+1)-1.           *      
!  INB    : vector INB(1:NNODES+1) with indices for the sublists *      
!           in NEIGHB. We must have that                         *      
!             INB(NNODES+1) = 1 + number of elements in NEIGHB . *      
!  IDEG   : vector IDEG(1:NNODES) of degrees of the nodes        *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  IROOT  : new starting node                                    *      
!  NLV    : number of levels                                     *      
!  LEVEL  : vector LEVEL(1:NNODES) listing nodes of identical    *      
!           level. For I=1, ..., NLV, the vector LEVEL contains  *      
!           the indices of the nodes of level I in positions     *      
!           LEVEL(K) for  K = ILV(I), ..., ILV(I+1)-1.           *      
!  ILV    : vector ILV(1:NLV+1) with indices for the level list  *      
!           LEVEL(). We must have  ILV(NLV+1) = 1 + number of    *      
!           entries in LEVEL .  If the graph is connected, this  *      
!           number is equal to  NNODES + 1.                      *      
!  LNODES : number of nodes in the component                     *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!  MARK   : LOGICAL MARK(1:NNODES), marking some nodes           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:                                         *      
!                                                                *      
!  LVSTRU    determines the level structure of one component of  *      
!            the graph                                           *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NEIGHB ( * ), INB ( * ), IDEG ( * ), LEVEL ( * ), ILV ( * &
      )                                                                 
      LOGICAL MARK ( * ) 
!                                                                       
      NLVOLD = 0 
   10 CONTINUE 
      CALL LVSTRU (IROOT, NNODES, NEIGHB, INB, MARK, NLV, LEVEL, ILV,   &
      LNODES)                                                           
!                                                                       
!        This is the exit from the loop                                 
!                                                                       
      IF (NLV.LE.NLVOLD) RETURN 
!                                                                       
      NLVOLD = NLV 
!                                                                       
!        Search for node of minimal degree in the previous level        
!                                                                       
      IMIN = ILV (NLV) 
      IDGMIN = IDEG (LEVEL (IMIN) ) 
      DO 20 I = ILV (NLV) + 1, ILV (NLV + 1) - 1 
         IF (IDEG (LEVEL (I) ) .LT.IDGMIN) THEN 
            IMIN = I 
            IDGMIN = IDEG (LEVEL (I) ) 
         ENDIF 
   20 END DO 
!                                                                       
!        Use this node as a start to construct                          
!        the level structure                                            
!                                                                       
      IROOT = LEVEL (IMIN) 
      GOTO 10 
      RETURN 
      END SUBROUTINE FNDROO                         
!                                                                       
!                                                                       
      INTEGER FUNCTION IBDWID (NNODES, NEIGHB, INB, NOLD, NNEW) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for the Cuthill-McKee algorithm.                   *      
!  It determines the band width of a matrix with known           *      
!  permutation of the nodes of its graph.                        *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NNODES : )  same as output from SUBROUTINE BLDGPH.            *      
!  NEIGHB : )  These vectors describe the graph of the matrix.   *      
!  INB    : )                                                    *      
!  NOLD   : permutation vector NOLD(1:NNODES)                    *      
!           For I=1, ..., NNODES,  NORG(I) designates the        *      
!           original index of the node now numbered by I.        *      
!  NNEW   : vector NNEW(1:NNODES), the inverse permutation of    *      
!           NOLD.                                                *      
!           For I=1, ..., NNODES,  NNEW(I) indicates the new     *      
!           index of the node previously numbered by I.          *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  IBDWID : Band width, i.e., the maximal distance of two        *      
!           adjacent nodes after the permutation                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NEIGHB ( * ), INB ( * ), NNEW ( * ), NOLD ( * ) 
!                                                                       
      MAXBW = 0 
      DO 10 I = 1, NNODES - 1 
         DO 20 K = INB (NOLD (I) ), INB (NOLD (I) + 1) - 1 
            IDIFF = ABS (I - NNEW (NEIGHB (K) ) ) 
            IF (IDIFF.GT.MAXBW) MAXBW = IDIFF 
   20    END DO 
   10 END DO 
      IBDWID = MAXBW 
      RETURN 
      END FUNCTION IBDWID                           
!                                                                       
!                                                                       
      SUBROUTINE ckbdwd (n, neighb, inb, m, icm, icmrev) 
!                                                                       
!*****************************************************************      
! check if by applying the Cuthill-McKee permutation the band    *      
! width might be enlarged instead of being reduced. In this      *      
! unlucky case the matrix shall remain unchanged, that is shall  *      
! not be permuted.                                               *      
!                                                                *      
! Input parameters:                                              *      
! =================                                              *      
! N          order of the sparse matrix          \  incidence    *      
! NEIGHB     indices of adjacent nodes            > graph        *      
! INB        indices for NEIGHB                  /  (see BLDGPH) *      
! M          half band width of the matrix produced by applying  *      
!            the Cuthill-McKee permutation                       *      
! ICM        (1:N)-vector with the CM-permutation. For i=1(1)n   *      
!            ICM(i) is the original number of the node with the  *      
!            new number i.                                       *      
! ICMREV     (1:N)-vector with the inverse permutation of ICM.   *      
!            For i=1(1)n ICMREV(i) is the new number of the node *      
!            with the old number i.                              *      
!                                                                *      
! Output parameters:                                             *      
! ==================                                             *      
! M          unchanged, if M <= original band width, else the    *      
!            original band width                                 *      
! ICM        (1:N)-vector with the identical permutation, if M   *      
!            had to be corrected, else with the CM-permutation   *      
! ICMREV     (1:N)-vector with the inverse permutation of ICM    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Juergen Dietel                                    *      
!  Date      : 06.19.1996                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER n, neighb ( * ), inb ( * ), m, icm ( * ), icmrev ( * ) 
!     half band width of the original matrix                            
      INTEGER morig 
      INTEGER diff 
      INTEGER i 
      INTEGER j 
!                                                                       
!                                                                       
      morig = 0 
      DO 20 i = 1, n 
         DO 10 j = inb (i), inb (i + 1) - 1 
            diff = abs (i - neighb (j) ) 
            IF (diff.gt.morig) morig = diff 
   10    END DO 
   20 END DO 
!                                                                       
!     new band width not smaller, but really bigger???                  
!     In this unlucky case the matrix shall remain unchanged.           
!                                                                       
      IF (m.gt.morig) then 
         m = morig 
         DO 30 i = 1, n 
            icm (i) = i 
            icmrev (i) = i 
   30    END DO 
      ENDIF 
                                                                        
      END SUBROUTINE ckbdwd                         
!                                                                       
!                                                                       
      SUBROUTINE LVSTRU (IROOT, NNODES, NEIGHB, INB, MARK, NLV, LEVEL,  &
      ILV, LNODES)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for the Cuthill-McKee algorithm.                   *      
!  It constructs the level structure of the component of a graph *      
!  generated by IROOT.                                           *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IROOT  : Number of the starting node                          *      
!  NNODES : Number of nodes of the graph. The graph is defined by*      
!           the following two vectors:                           *      
!  NEIGHB : vector NEIGHB(1:*) with the lists of adjacent nodes  *      
!           For I=1, ..., NNODES, the vector NEIGHB contains the *      
!           indices of the nodes that are adjacent to node I in  *      
!           NEIGHB(K) for  K=INB(I), ...,  INB(I+1)-1.           *      
!  INB    : vector INB(1:NNODES+1) with indices for the sublists *      
!           in NEIGHB. We must have that                         *      
!             INB(NNODES+1) = 1 + number of elements in NEIGHB . *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NLV    : number of levels                                     *      
!  LEVEL  : vector LEVEL(1:NNODES) listing nodes of identical    *      
!           level. For I=1, ..., NLV, the vector LEVEL contains  *      
!           the indices of the nodes of level I in positions     *      
!           LEVEL(K) for  K = ILV(I), ..., ILV(I+1)-1.           *      
!  ILV    : vector ILV(1:NLV+1) with indices for the level list  *      
!           LEVEL(). We must have  ILV(NLV+1) = 1 + number of    *      
!           entries in LEVEL .  If the graph is connected, this  *      
!           number is equal to  NNODES + 1.                      *      
!  LNODES : number of nodes in the component                     *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!  MARK   : LOGICAL MARK(1:NNODES), marking some nodes           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NEIGHB ( * ), INB (1:NNODES + 1), LEVEL (1:NNODES),       &
      ILV ( * )                                                         
      LOGICAL MARK (1:NNODES) 
!                                                                       
      NLV = 0 
      LEVEL (1) = IROOT 
      MARK (IROOT) = .TRUE. 
      N = 1 
!                                                                       
!     N is the number of nodes already found                            
!                                                                       
      LEVEND = 0 
   10 CONTINUE 
      NLV = NLV + 1 
!                                                                       
!        NLV is the number of already found levels                      
!                                                                       
      LEVBEG = LEVEND+1 
      LEVEND = N 
!                                                                       
!        LEVBEG denotes the start, while LEVEND denotes the end         
!        of the last level in LEVEL()                                   
!                                                                       
      ILV (NLV) = LEVBEG 
!                                                                       
!        Find nodes of level  NLV+1:                                    
!        search for nodes adjacent to nodes of level NLV,               
!        record them in LEVEL, if they have not been marked             
!                                                                       
      DO 20 I = LEVBEG, LEVEND 
         DO 30 J = INB (LEVEL (I) ), INB (LEVEL (I) + 1) - 1 
            IF (.NOT.MARK (NEIGHB (J) ) ) THEN 
               N = N + 1 
               LEVEL (N) = NEIGHB (J) 
               MARK (NEIGHB (J) ) = .TRUE. 
            ENDIF 
   30    END DO 
   20 END DO 
      IF (N.GT.LEVEND) GOTO 10 
!                                                                       
      LNODES = LEVEND 
      ILV (NLV + 1) = LNODES + 1 
!                                                                       
!     Mark all recently marked nodes .FALSE.                            
!                                                                       
      DO 40 I = 1, LNODES 
         MARK (LEVEL (I) ) = .FALSE. 
   40 END DO 
      RETURN 
      END SUBROUTINE LVSTRU                         
!                                                                       
!                                                                       
      SUBROUTINE PERMUT (N, IPERM, XOLD, XNEW) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for solving linear systems of equations with the   *      
!  Cuthill-McKee algorithm.                                      *      
!  PERMUT performs the permutation defined in IPERM on XOLD.     *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : order of the vector                                  *      
!  IPERM  : vector IPERM(1:N) with the permutations of 1, ..., N *      
!  XOLD   : vector XOLD(1:N) whose entries are to be permuted    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  XNEW   : vector XNEW(1:N), the permuted vector                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 17.11.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      INTEGER IPERM ( * ) 
      DIMENSION XOLD ( * ), XNEW ( * ) 
!                                                                       
      DO 10 I = 1, N 
         XNEW (IPERM (I) ) = XOLD (I) 
   10 END DO 
      RETURN 
      END SUBROUTINE PERMUT                         
!                                                                       
!                                                                       
      SUBROUTINE RDMTRX (FLNAM, MAXELM, MAXROW, NROW, NV, V, IC, IR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for solving of linear systems with sparse matrices.*      
!  It opens a data file and reads matrix elements into the vector*      
!  V. It thus prepares for the Cuthill-McKee algorithm.          *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FLNAM  : Name of the input file with the following structure: *      
!           - number of rows of the matrix                       *      
!           - for each row:                                      *      
!             - for each nonzero element we store the tupel:     *      
!               <column index> <element>                         *      
!               (indices start at 1)                             *      
!             - a tupel with column index 0 denotes the end of   *      
!               the row                                          *      
!  MAXELM : Dimensioning number for arrays, that will contain    *      
!           matrix elements. MAXELM must be at least equal to the*      
!           number of nonzero matrix elements.                   *      
!  MAXROW : Dimensioning number for auxiliary arrays.            *      
!           MAXROW must be at least as large as the number of    *      
!           matrix rows.                                         *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NROW   : number of rows                                       *      
!  NV     : number of entries in the vector V ( = number of non- *      
!           zero entries in the original matrix). The vector IC  *      
!           has the same number of entries.                      *      
!  V      : vector V(1:NV) with nonzero entries. The nonzero     *      
!           elements of the original matrix are written row-wise *      
!           and in sequence in V.                                *      
!  IC     : vector IC(1:NV) of column indices. For each V(I), the*      
!           value in IC(I) denotes the column of the matrix that *      
!           contained  V(I). (Indices start with 1)              *      
!  IR     : vector IR(1:NROW+1) with row pointers. For I=1, ..., *      
!           , NROW,  IR(I) denotes the index, where row  I starts*      
!           in  V and IC. We must have  IR(NROW+1) =  NV+1.      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
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
      INTEGER IC ( * ), IR ( * ) 
      DIMENSION V ( * ) 
      CHARACTER ( * ) FLNAM 
!                                                                       
      PARAMETER (IFLNO = 8) 
!                                                                       
      OPEN (UNIT = IFLNO, FILE = FLNAM) 
      READ (IFLNO, * ) NROW 
      IF (nrow + 1.gt.maxrow) then 
         WRITE ( * , * ) 'RDMTRX: Matrix too big (maximum: ', maxrow -  &
         1, 'rows)!'                                                    
         STOP 
      ENDIF 
      NV = 0 
      DO 10 I = 1, NROW 
         IR (I) = NV + 1 
   20    CONTINUE 
         READ (IFLNO, * ) IC0, V0 
         IF (IC0.EQ.0) GOTO 10 
         NV = NV + 1 
         IF (nv.gt.maxelm) then 
            WRITE ( * , * ) 'RDMTRX: Matrix has too many non-zero',     &
            'elements (maximum: ', maxelm, ')!'                         
            STOP 
         ENDIF 
         V (NV) = V0 
         IC (NV) = IC0 
         GOTO 20 
   10 END DO 
      CLOSE (IFLNO) 
      IR (NROW + 1) = NV + 1 
      RETURN 
      END SUBROUTINE RDMTRX                         
!                                                                       
!                                                                       
      SUBROUTINE SRTDEG (NODE, IDEG, IBEG, IEND) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine for  SUBROUTINE CUTHIL.                            *      
!  Sorts a partial list of nodes in the vector NODE according to *      
!  increasing degrees.                                           *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NODE   : vector NODE(1:*) with node numbers. Each entry must  *      
!           be an index of IDEG.                                 *      
!  IDEG   : vector IDEG(1:*) with the node degrees. The node     *      
!           NODE(I) has the degree   IDEG(NODE(I)).              *      
!  IBEG   : starting index for the nodes that shall be sorted.   *      
!  IEND   : final index of the nodes to be sorted                *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  NODE   : as on input, except for partial reordering           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.17.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NODE ( * ), IDEG ( * ) 
!                                                                       
      DO 10 I = IBEG + 1, IEND 
         NODE0 = NODE (I) 
         IDEG0 = IDEG (NODE0) 
!                                                                       
         K0 = I 
         DO 20 K = I - 1, IBEG, - 1 
            IF (IDEG0.GE.IDEG (NODE (K) ) ) GOTO 30 
            NODE (K + 1) = NODE (K) 
            K0 = K 
   20    END DO 
   30    NODE (K0) = NODE0 
   10 END DO 
      RETURN 
      END SUBROUTINE SRTDEG                         

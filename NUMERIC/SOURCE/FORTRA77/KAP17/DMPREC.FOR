C[BA*)
C[LE*)
      SUBROUTINE DMPREC (DSMALL,DLARGE)
C[IX{DMPREC}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the machine constant in DOUBLE      *
C PRECISION and determines the largest representable DOUBLE      *
C PRECISION number for the machine used.                         *
C See the description of DLARGE !                                *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETER: none                                          *
C ================                                               *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C DSMALL  : smallest DOUBLE PRECISION number for which           *
C              1.D0 + DSMALL > 1.D0    holds.                *
C DLARGE  : largest representable DOUBLE PRECISION number.       *
C           DLARGE is given as a constant. The assignment        *
C           DLARGE = .... must be an executable FORTRAN          *
C           statement. The other given DLARGE assignment must be *
C           commented out. If DLARGE is not given for the used   *
C           computer, then its proper value must be computed and *
C           the programm must be amended accordingly.            *
C           As an example we initialize DLARGE for the CYBER 930.*
C                                                                *
C                                                                *
C LOCKAL VARIABLES: none                                         *
C ================                                               *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker KrÅger                                      *
C  Date     : 07.08.1990                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declare DSMALL, DLARGE DOUBLE PRECISION
C
      DOUBLE PRECISION DSMALL,DLARGE
C
C Initialize DSMALL
C
      DSMALL=1.0D0
C
C Compute DSMALL
C      While  1 < 1 + DSMALL  holds,
C      we half DSMALL.
C
10    DSMALL=5.0D-01*DSMALL
      IF(MACHPD(1.0D0+DSMALL) .EQ. 1) GOTO 10
C
C  DSMALL is the smallest DOUBLE PRECISION number with
C  1 + DSMALL > 1
C
      DSMALL=2.0D0*DSMALL
C
C For the CYBER 930 (for example)
C
C      DLARGE=3.03957169355461521391917692D+1232
C
C For the PC/AT MS-FORTRAN  (for example)
C
      DLARGE=1.046395124205339D+308
C
      RETURN
      END

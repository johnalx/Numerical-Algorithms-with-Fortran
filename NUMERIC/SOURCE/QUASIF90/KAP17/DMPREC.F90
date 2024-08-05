      SUBROUTINE DMPREC (DSMALL, DLARGE) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the machine constant in DOUBLE      *      
! PRECISION and determines the largest representable DOUBLE      *      
! PRECISION number for the machine used.                         *      
! See the description of DLARGE !                                *      
!                                                                *      
!                                                                *      
! INPUT PARAMETER: none                                          *      
! ================                                               *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! DSMALL  : smallest DOUBLE PRECISION number for which           *      
!              1.D0 + DSMALL > 1.D0    holds.                *          
! DLARGE  : largest representable DOUBLE PRECISION number.       *      
!           DLARGE is given as a constant. The assignment        *      
!           DLARGE = .... must be an executable FORTRAN          *      
!           statement. The other given DLARGE assignment must be *      
!           commented out. If DLARGE is not given for the used   *      
!           computer, then its proper value must be computed and *      
!           the programm must be amended accordingly.            *      
!           As an example we initialize DLARGE for the CYBER 930.*      
!                                                                *      
!                                                                *      
! LOCKAL VARIABLES: none                                         *      
! ================                                               *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker KrÅger                                      *      
!  Date     : 07.08.1990                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declare DSMALL, DLARGE DOUBLE PRECISION                               
!                                                                       
      DOUBLEPRECISION DSMALL, DLARGE 
!                                                                       
! Initialize DSMALL                                                     
!                                                                       
      DSMALL = 1.0D0 
!                                                                       
! Compute DSMALL                                                        
!      While  1 < 1 + DSMALL  holds,                                    
!      we half DSMALL.                                                  
!                                                                       
   10 DSMALL = 5.0D-01 * DSMALL 
      IF (MACHPD (1.0D0 + DSMALL) .EQ.1) GOTO 10 
!                                                                       
!  DSMALL is the smallest DOUBLE PRECISION number with                  
!  1 + DSMALL > 1                                                       
!                                                                       
      DSMALL = 2.0D0 * DSMALL 
!                                                                       
! For the CYBER 930 (for example)                                       
!                                                                       
!      DLARGE=3.03957169355461521391917692D+1232                        
!                                                                       
! For the PC/AT MS-FORTRAN  (for example)                               
!                                                                       
      DLARGE = 1.046395124205339D+308 
!                                                                       
      RETURN 
      END SUBROUTINE DMPREC                         

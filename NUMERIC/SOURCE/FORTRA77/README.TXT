additional files for working with GNU FORTRAN 0.5.18 under OS/2:

makefile.mk     default makefile for DMAKE 3.80
makefile.in     input file for call of DMAKE as follows:
                `dmake test <makefile.in >makefile.aus'
os2dmake.cmd    call of DMAKE keeping the compiler in memory for 30 min
                in order to increase speed of compilation
gnutst.cmd      complete test of EFNUM including call of DMAKE and
                redirection of any output into a file

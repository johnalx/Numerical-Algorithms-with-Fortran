# makefile for EFNUM 8.04 and DMAKE 3.80
# Author: Juergen Dietel, Computer Center, RWTH Aachen
# Call:    dmake                -    like `dmake all'
#          dmake all            -    create all test programs
#          dmake mpgauss.exe    -    create only `mpgauss.exe'
#          ...                  -    ...
#          dmake test           -    run all test programs
#          dmake clean          -    delete all created files
# Add-ons: DEB=j                -    include debugging information
#          OPT=j                -    turn on compiler optimization
#          VERB=j               -    compile and link verbosely
#          EMX09BOLD=j          -    emx 0.9b with GNU FORTRAN 0.5.17
#                                    (default: emx 0.9b with
#                                    GNU FORTRAN 0.5.18 or maybe even
#                                    the older combination of emx 0.9a
#                                    and GNU FORTRAN 0.5.17)

# macros for compiling EFNUM with GNU FORTRAN 0.5.18 under OS/2
FC     =          g77
OUTPUT =          -o $@
CFLAGS =          -c $(OUTPUT) $(DEBUG) $(OPTIM) $(VERBOS) #-Wall
LFLAGS =          $(OUTPUT)
LINK   =          $(LINKER) $(LFLAGS) $& $(LFLAG2)

.IF $(DEB) == j
DEBUG  =          -g
.ENDIF

.IF $(OPT) == j
OPTIM  =          -O
.ENDIF

.IF $(VERB) == j
VERBOS =          -v
.ENDIF

.IF $(EMX09BOLD) == $(NULL)
# needed for emx 0.9a with GNU CC 2.6.3 and GNU FORTRAN 0.5.17 or
# for emx 0.9b with GNU CC 2.7.2 and GNU FORTRAN 0.5.18:
LINKER =          $(FC)
.ELSE
# needed for emx 0.9b with GNU CC 2.7.2 and GNU FORTRAN 0.5.17, because
# without the library file `c_alias.a' the linker reports errors. The
# call `g77' does not work, as it causes the libraries to be searched in
# wrong order. So you have to go back to `gcc' and explicitly name all
# libraries needed additionally for FORTRAN.
LINKER =          gcc
LFLAG2 =          -lf2c -lm -lc_alias
.ENDIF

# macro names for all folders containing EFNUM source files
00     =          kap00
02     =          kap02
03     =          kap03
04     =          kap04
05     =          kap05
06     =          kap06
07     =          kap07
08     =          kap08
09     =          kap09
10     =          kap10
11     =          kap11
12     =          kap12
13     =          kap13
14     =          kap14
15     =          kap15
16     =          kap16
17     =          kap17
18     =          kap18
02t    =          $(02)\tstpro
03t    =          $(03)\tstpro
04t    =          $(04)\tstpro
05t    =          $(05)\tstpro
06t    =          $(06)\tstpro
07t    =          $(07)\tstpro
08t    =          $(08)\tstpro
09t    =          $(09)\tstpro
10t    =          $(10)\tstpro
11t    =          $(11)\tstpro
12t    =          $(12)\tstpro
13t    =          $(13)\tstpro
14t    =          $(14)\tstpro
15t    =          $(15)\tstpro
16t    =          $(16)\tstpro
17t    =          $(17)\tstpro
18t    =          $(18)\tstpro

# lists with EFNUM file names (without extension `.for')
# all files with methods from chapter 4
l04    =          gauss gausrs choky chonpd cg pivot trdig trdsy trdnpd\
                  cyctr cytsy cytnpd fdiag fdisy band chobnd shouse\
                  hacond cline condae postit cutcho cutgau
# all files with test programs for chapter 4
l04t   =          mpgauss mpgausrs mpchoky mpchonpd mpcg mppivot\
                  mptrdig mptrdsy mptrdnpd mpcyctr mpcytnpd mpfdiag\
                  mpfdisy mpband mpchobnd mpshouse mpcline mpcondae\
                  mpcutcho mpcutgau
# all files with methods from chapter 8
l08    =          gadesm polfit slfit slpre senorm rfft fft fftb fdicht\
                  fourn ffako ffakon snlfit snlpre
# all files with test programs for chapter 8
l08t   =          mpgadesm mppolfit mpslfit mprfft mpfft mpfftb\
                  mpfdicht mpfourn mpffako mpffakon mpsnlfit
# all files with methods from chapter 10
l10    =          isplnp ispl1d ispl2d ispl3d isplpe isplnk sptab spval\
                  isplpa psppv psptab pspval ispltr tsptab tspval\
                  hermit hmtab parmit pmtab
# all files with test programs for chapter 10
l10t   =          misplnp1 misplnp2 misplnp3 misplpa1 misplpa2 mpispltr\
                  mphermit mpparmit mptspana mppspana
# all files with methods from chapter 11
l11    =          cfspnp cfsp1d cfsp2d cfsp3d cfsppe ncyfsy cfsppa\
                  cfsptr
# all files with methods from chapter 12
l12    =          bicsp1 bicsp2 bicsp3 bspval prob2 pr2tra prob3 pr3tra\
                  sesspm cepspm zspmok cubbez bezier calcp calcvp\
                  calcwp mocube curvp knotvo deboor curvpg knotvg bsplfl
# all files with test programs for chapter 12
l12t   =          mpbicsp1 mpbicsp2 mpbicsp3 mpprob2 mppr2tra mpprob3\
                  mppr3tra mpbezier mp2bezie mpmocube mpcurvp mpcurvpg\
                  mpbsplfl
# all files with methods from chapter 15
l15    =          quanec ortogp ccferr clensh wgknot quarom gax gale0
# all files with test programs for chapter 15
l15t   =          mpquanec mportogp mpccferr mpquarom mpgax
# all files with methods from chapter 16
l16    =          k4necn k3nec3 k4buri k4rori k3rori k4gaue k4gauv\
                  k3gaun fibicu fibic2
# all files with test programs for chapter 16
l16t   =          mpk4necn mpk3nec3 mpk4buri mpk4rori mpk3rori mpk4gaue\
                  mpk4gauv mpk3gaun mpfibicu mpfibic2
# alle prerequisite files for `mprktrb.exe' from $(l17)
l17_0  =          rktrb dmprec\
                  coeffi hstart hull civps rkf43 rkf54 rk546m rke54\
                  hiha5 rk547m rk547s rk547c rk658m rk658s rk658c rkv65\
                  rkf65a rkf65b rkc65 rkv65a rkv65b rkv76 rk8713 rkf87\
                  rkv87 rkv98 vmnorm rkstep step32
# all prerequisite files for `mpbvp.exe' from $(l17)
l17_1  =          ivp desabm desext irkdrv irkcoe imruku $(l17_0)
# all files with methods from chapter 17
l17    =          deqssp frkfsy gear4 $(l17_1)
# all files with test programs for chapter 17
l17t   =          mpdeqssp mprktrb mpirkdrv mpivp mpfrkfsy mpdesabm\
                  mpdesext mpgear4
# all prerequisite files for `mprktrb.exe' (excluding `mprktrb.obj')
OBJS0  =          machpd $(l17_0)
# all prerequisite files for `mpbvp.exe' (excluding `mpbvp.obj')
OBJS1  =          bvp gauss machpd gale0 $(l17_1)

# liste of all executable files (without extension `.exe' yet)
EXES   =          mpnewton mppegasu mpzerorf mpzeroin\
                  mpmullrp mpbaupol mplaguer\
                  mpgauss mpgausrs mpchoky mpchonpd mpcg mppivot\
                  mptrdig mptrdsy mptrdnpd mpcyctr mpcytnpd mpfdiag\
                  mpfdisy mpband mpchobnd mpshouse mpcline mpcondae\
                  mpcutcho mpcutgau\
                  mpadsor\
                  mpsmnewd mpsmnewt mpnlesys mpbrown\
                  mpeigen mpeigval\
                  mpgadesm mppolfit mpslfit mprfft mpfft\
                  mpfftb mpfdicht mpfourn mpffako mpffakon mpsnlfit\
                  mpnewtip mpralgip mpshpglo mpshplok mpshpflg\
                  misplnp1 misplnp2 misplnp3 misplpa1 misplpa2 mpispltr\
                  mphermit mpparmit mptspana mppspana\
                  mpcfspnp mpcfpab mpcfpap mpcfsppa mpcfsptr\
                  mpbicsp1 mpbicsp2 mpbicsp3 mpprob2 mppr2tra mpprob3\
                  mppr3tra mpbezier mp2bezie mpmocube mpcurvp mpcurvpg\
                  mpbsplfl\
                  mpakima1 mpakima2 mprenn2d mprenn3d\
                  mpdifrom\
                  mpquanec mportogp mpccferr mpquarom mpgax\
                  mpk4necn mpk3nec3 mpk4buri mpk4rori mpk3rori mpk4gaue\
                  mpk4gauv mpk3gaun mpfibicu mpfibic2\
                  mpdeqssp mprktrb mpirkdrv mpivp mpfrkfsy mpdesabm\
                  mpdesext mpgear4\
                  mpbvp

# association of EFNUM source file names with folders containing the
# files
.SETDIR=$(00):    {machpd}.for
.SETDIR=$(02):    {newpsz newpol newmod pegasu zerorf zeroin}.for
.SETDIR=$(02t):   {mpnewton mppegasu mpzerorf mpzeroin}.for
.SETDIR=$(03):    {muller baupol laguer}.for
.SETDIR=$(03t):   {mpmullrp mpbaupol mplaguer}.for
.SETDIR=$(04):    {$(l04)}.for
.SETDIR=$(04t):   {$(l04t)}.for
.SETDIR=$(05):    {adsor}.for
.SETDIR=$(05t):   {mpadsor}.for
.SETDIR=$(06):    {smnewd smnewt fenorm nlesys jacobi brown}.for
.SETDIR=$(06t):   {mpsmnewd mpsmnewt mpnlesys mpbrown}.for
.SETDIR=$(07):    {eigval eigen}.for
.SETDIR=$(07t):   {mpeigen mpeigval}.for
.SETDIR=$(08):    {$(l08)}.for
.SETDIR=$(08t):   {$(l08t)}.for
.SETDIR=$(09):    {newtip ralgip shpglo shplok shpflg}.for
.SETDIR=$(09t):   {mpnewtip mpralgip mpshpglo mpshplok mpshpflg}.for
.SETDIR=$(10):    {$(l10)}.for
.SETDIR=$(10t):   {$(l10t)}.for
.SETDIR=$(11):    {$(l11)}.for
.SETDIR=$(11t):   {mpcfspnp mpcfpab mpcfpap mpcfsppa mpcfsptr}.for
.SETDIR=$(12):    {$(l12)}.for
.SETDIR=$(12t):   {$(l12t)}.for
.SETDIR=$(13):    {akima1 akima2 renn2d renn3d}.for
.SETDIR=$(13t):   {mpakima1 mpakima2 mprenn2d mprenn3d}.for
.SETDIR=$(14):    {splfvd tspana pspana difrom}.for
.SETDIR=$(14t):   {mpdifrom}.for
.SETDIR=$(15):    {$(l15)}.for
.SETDIR=$(15t):   {$(l15t)}.for
.SETDIR=$(16):    {$(l16)}.for
.SETDIR=$(16t):   {$(l16t)}.for
.SETDIR=$(17):    {$(l17)}.for
.SETDIR=$(17t):   {$(l17t)}.for
.SETDIR=$(18):    {bvp}.for
.SETDIR=$(18t):   {mpbvp}.for

# implicit rule: OBJ file --> corresponding EXE file
%.exe:            %.obj
	%$(LINK)

# implicit rule: FOR file --> corresponding OBJ file
%.obj .PRECIOUS:  %.for
	$(FC) $(CFLAGS) $<

# pseudo target `all' for all test programms together (default target)
all:              {$(EXES)}.exe

# pseudo target `test' for calling all test programs
test:             all
	mpnewton
	mppegasu
	mpzerorf
	mpzeroin
	mpmullrp
	mpbaupol
	mplaguer
	mpgauss
	mpgausrs
	mpchoky
	mpchonpd
	mpcg
	mppivot
	mptrdig
	mptrdsy
	mptrdnpd
	mpcyctr
	mpcytnpd
	mpfdiag
	mpfdisy
	mpband
	mpchobnd
	mpshouse
	mpcline
	mpcondae
	+(cd kap04\tstpro & ..\..\mpcutcho & ..\..\mpcutgau & cd ..\..)
	mpadsor
	mpsmnewd
	mpsmnewt
	mpnlesys
	mpbrown
	mpeigen
	mpeigval
	mpgadesm
	mppolfit
	mpslfit
	mprfft
	mpfft
	mpfftb
	mpfdicht
	mpfourn
	mpffako
	mpffakon
	mpsnlfit
	mpnewtip
	mpralgip
	mpshpglo
	mpshplok
	mpshpflg
	misplnp1
	misplnp2
	misplnp3
	misplpa1
	misplpa2
	mpispltr
	mphermit
	mpparmit
	mptspana
	mppspana
	mpcfspnp
	mpcfpab
	mpcfpap
	mpcfsppa
	mpcfsptr
	mpbicsp1
	mpbicsp2
	mpbicsp3
	mpprob2
	mppr2tra
	mpprob3
	mppr3tra
	mpbezier
	+(cd kap12\tstpro & ..\..\mp2bezie & cd ..\..)
	mpmocube
	mpcurvp
	mpcurvpg
	mpbsplfl
	mpakima1
	mpakima2
	mprenn2d
	mprenn3d
	mpdifrom
	mpquanec
	mportogp
	mpccferr
	mpquarom
	mpgax
	mpk4necn
	mpk3nec3
	mpk4buri
	mpk4rori
	mpk3rori
	mpk4gaue
	mpk4gauv
	mpk3gaun
	mpfibicu
	mpfibic2
	mpdeqssp
	mprktrb
	mpirkdrv
	mpivp
	mpfrkfsy
	mpdesabm
	mpdesext
	mpgear4
	mpbvp

# additional rules for all test programs from chapter 2
mpnewton.exe:     {newpol newpsz machpd newmod}.obj
mppegasu.exe:     {pegasu machpd}.obj
mpzerorf.exe:     {zerorf machpd}.obj
mpzeroin.exe:     {zeroin machpd}.obj

# additional rules for all test programs from chapter 3
mpmullrp.exe:     {muller}.obj
mpbaupol.exe:     {baupol}.obj
mplaguer.exe:     {laguer machpd}.obj

# additional rules for all test programs from chapter 4
mpgauss.exe:      {gauss machpd hacond postit}.obj
mpgausrs.exe:     {gausrs gauss machpd}.obj
mpchoky.exe:      {choky machpd}.obj
mpchonpd.exe:     {chonpd choky machpd}.obj
mpcg.exe:         {cg machpd}.obj
mppivot.exe:      {pivot machpd}.obj
mptrdig.exe:      {trdig machpd}.obj
mptrdsy.exe:      {trdsy machpd}.obj
mptrdnpd.exe:     {trdnpd trdsy machpd}.obj
mpcyctr.exe:      {cyctr machpd}.obj
mpcytnpd.exe:     {cytnpd cytsy machpd}.obj
mpfdiag.exe:      {fdiag machpd}.obj
mpfdisy.exe:      {fdisy machpd}.obj
mpband.exe:       {band machpd}.obj
mpchobnd.exe:     {chobnd machpd}.obj
mpshouse.exe:     {shouse machpd}.obj
mpcline.exe:      {cline gauss machpd}.obj
mpcondae.exe:     {condae machpd gauss}.obj
mpcutcho.exe:     {cutcho chobnd machpd}.obj
mpcutgau.exe:     {cutgau cutcho band chobnd machpd}.obj

# additional rules for all test programs from chapter 5
mpadsor.exe:      {adsor machpd}.obj

# additional rules for all test programs from chapter 6
mpsmnewd.exe:     {smnewd fenorm gauss machpd}.obj
mpsmnewt.exe:     {smnewt fenorm gauss machpd}.obj
mpnlesys.exe:     {nlesys machpd senorm jacobi shouse}.obj
mpbrown.exe:      {brown jacobi}.obj

# additional rules for all test programs from chapter 7
mpeigval.exe:     {eigval}.obj
mpeigen.exe:      {eigen}.obj

# additional rules for all test programs from chapter 8
mpgadesm.exe:     {gadesm choky machpd}.obj
mppolfit.exe:     polfit.obj
mpslfit.exe:      {slfit slpre shouse senorm machpd}.obj
mprfft.exe:       {rfft}.obj
mpfft.exe:        {fft}.obj
mpfftb.exe:       {fft fftb}.obj
mpfdicht.exe:     {fdicht fft fftb}.obj
mpfourn.exe:      {fourn fft fftb}.obj
mpffako.exe:      {ffako fft fftb}.obj
mpffakon.exe:     {ffakon fft fftb}.obj
mpsnlfit.exe:     {snlfit machpd senorm snlpre shouse}.obj

# additional rules for all test programs from chapter 9
mpnewtip.exe:     {newtip}.obj
mpralgip.exe:     {ralgip}.obj
mpshpglo.exe:     {shpglo}.obj
mpshplok.exe:     {shplok}.obj
mpshpflg.exe:     {shpflg}.obj

# additional rules for all test programs from chapter 10
misplnp1.exe:     {ispl1d trdsy machpd isplnp ispl2d ispl3d}.obj \
                  {isplpe isplnk cytsy trdig sptab}.obj
misplnp2.exe:     {isplnp ispl1d ispl2d ispl3d isplpe isplnk}.obj \
                  {trdsy cytsy trdig machpd sptab}.obj
misplnp3.exe:     {isplnp ispl1d ispl2d ispl3d isplpe isplnk}.obj \
                  {trdsy cytsy trdig machpd sptab}.obj
misplpa1.exe:     {isplpa psppv ispl1d ispl2d isplpe trdsy cytsy}.obj \
                  {machpd psptab}.obj
misplpa2.exe:     {isplpa psppv ispl1d ispl2d isplpe trdsy cytsy}.obj \
                  {machpd psptab}.obj
mpispltr.exe:     {ispltr isplpe cytsy machpd tsptab}.obj
mphermit.exe:     {hermit cyctr machpd hmtab}.obj
mpparmit.exe:     {parmit hermit cyctr machpd pmtab}.obj
mptspana.exe:     {ispltr isplpe cytsy machpd tspana splfvd}.obj
mppspana.exe:     {isplpa psppv ispl1d ispl2d isplpe trdsy cytsy}.obj \
                  {machpd pspana splfvd}.obj

# additional rules for all test programs from chapter 11
mpcfspnp.exe:     {cfspnp cfsp1d cfsp2d cfsp3d cfsppe fdisy fdiag}.obj \
                  {ncyfsy machpd sptab}.obj
mpcfpab.exe:      {cfsppa psppv cfsp1d cfsp2d cfsppe fdisy ncyfsy}.obj \
                  {machpd psptab}.obj
mpcfpap.exe:      {cfsppa psppv cfsp1d cfsp2d cfsppe fdisy ncyfsy}.obj \
                  {machpd psptab}.obj
mpcfsppa.exe:     {cfsppa psppv cfsp1d cfsp2d cfsppe fdisy ncyfsy}.obj \
                  {machpd cfspnp cfsp3d fdiag sptab}.obj
mpcfsptr.exe:     {cfsptr cfsppe ncyfsy machpd tsptab}.obj

# additional rules for all test programs from chapter 12
mpbicsp1.exe:     {bicsp1 hermit cyctr machpd bspval}.obj
mpbicsp2.exe:     {bicsp2 bicsp1 hermit cyctr machpd bspval}.obj
mpbicsp3.exe:     {bicsp3 bicsp2 bicsp1 hermit cyctr machpd bspval}.obj
mpprob2.exe:      {prob2 cepspm sesspm}.obj
mppr2tra.exe:     {pr2tra prob2 cepspm sesspm}.obj
mpprob3.exe:      {prob3 cepspm sesspm}.obj
mppr3tra.exe:     {pr3tra prob3 cepspm sesspm}.obj
mpbezier.exe:     {bezier calcvp calcp calcwp}.obj
mp2bezie.exe:     {bezier calcvp calcp}.obj
mpmocube.exe:     {mocube}.obj
mpcurvp.exe:      {curvp deboor knotvo}.obj
mpcurvpg.exe:     {curvpg deboor knotvg}.obj
mpbsplfl.exe:     {bsplfl deboor knotvo}.obj

# additional rules for all test programs from chapter 13
mpakima1.exe:     {akima1 machpd}.obj
mpakima2.exe:     {akima2 machpd}.obj
mprenn2d.exe:     {renn2d machpd}.obj
mprenn3d.exe:     {renn3d machpd}.obj

# additional rules for all test programs from chapter 14
mpdifrom.exe:     {difrom machpd}.obj

# additional rules for all test programs from chapter 15
mpquanec.exe:     {quanec}.obj
mportogp.exe:     {ortogp gauss muller machpd}.obj
mpccferr.exe:     {ccferr wgknot clensh}.obj
mpquarom.exe:     {quarom}.obj
mpgax.exe:        {gax gale0 wgknot clensh quarom quanec machpd}.obj

# additional rules for all test programs from chapter 16
mpk4necn.exe:     {k4necn}.obj
mpk3nec3.exe:     {k3nec3 k4rori}.obj
mpk4buri.exe:     {k4buri}.obj
mpk4rori.exe:     {k4rori}.obj
mpk3rori.exe:     {k3rori k3nec3 k4rori}.obj
mpk4gaue.exe:     {k4gaue}.obj
mpk4gauv.exe:     {k4gauv k4gaue}.obj
mpk3gaun.exe:     {k3gaun}.obj
mpfibicu.exe:     {fibicu bicsp2 bicsp1 hermit cyctr machpd}.obj
mpfibic2.exe:     {fibic2 bspval fibicu bicsp3 bicsp2 bicsp1}.obj \
                  {hermit cyctr machpd}.obj

# additional rules for all test programs from chapter 17
mpdeqssp.exe:     {deqssp machpd}.obj
mprktrb.exe:      {$(OBJS0)}.obj
mpirkdrv.exe:     {irkdrv irkcoe machpd imruku gale0}.obj
mpivp.exe:        {ivp machpd}.obj
mpfrkfsy.exe:     {frkfsy machpd}.obj
mpdesabm.exe:     {desabm machpd}.obj
mpdesext.exe:     {desext machpd}.obj
mpgear4.exe:      {gear4 machpd ivp gauss}.obj

# additional rules for all test programs from chapter 18
mpbvp.exe:        {$(OBJS1)}.obj

# pseudo target `clean' for tidying up
clean:
	- del *.obj
	- del *.exe
	- del *.au*

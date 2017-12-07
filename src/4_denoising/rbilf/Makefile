CFLAGS   = -std=c99 -O3 -DNDEBUG -w
CXXFLAGS = -O3 -DNDEBUG -w
CPPFLAGS = -Ilib/iio -Ilib/argparse -Ilib/tvl1flow
LDLIBS   = -ljpeg -lpng -ltiff -lm
OMPFLAGS = -fopenmp

BINDIR   = build/bin/
RNLM     = $(BINDIR)vnlmeans
TVL1     = $(BINDIR)tvl1flow
SCRIPT   = $(BINDIR)vnlm-gt.sh $(BINDIR)vnlm-mp4.sh $(BINDIR)psnr.sh
OBJ_RNLM = lib/iio/iio.o lib/argparse/argparse.o src/main.o
OBJ_TVL1 = lib/iio/iio.o lib/tvl1flow/main.o

all          : $(RNLM) $(TVL1) $(SCRIPT)
$(RNLM)      : $(OBJ_RNLM) $(BINDIR) ; $(CC) $(LDFLAGS) -o $@ $(OBJ_RNLM) $(LDLIBS)
$(TVL1)      : $(OBJ_TVL1) $(BINDIR) ; $(CC) $(LDFLAGS) $(OMPFLAGS) -o $@ $(OBJ_TVL1) $(LDLIBS)
$(BINDIR)%.sh: scripts/%.sh $(BINDIR); cp $< $@ 
$(BINDIR)    : ; mkdir -p $@
clean        : ; $(RM) $(NLDCT) $(TVL1) $(OBJ_TVL1) $(OBJ_RNLM) $(SCRIPT)


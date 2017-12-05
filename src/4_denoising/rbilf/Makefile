CFLAGS   = -std=c99 -O3 -DNDEBUG -w
CXXFLAGS = -O3 -DNDEBUG -w
CPPFLAGS = -Ilib/iio -Ilib/argparse
LDLIBS   = -ljpeg -lpng -ltiff -lm

BINDIR   = build/bin/
NLDCT    = $(BINDIR)vnlmeans
TVL1     = $(BINDIR)tvl1flow
SCRIPT   = $(BINDIR)vnlm-gt.sh $(BINDIR)vnlm-mp4.sh
OBJ      = lib/iio/iio.o lib/argparse/argparse.o src/main.o

all          : $(NLDCT) $(TVL1) $(SCRIPT)
$(NLDCT)     : $(OBJ) $(BINDIR) ; $(CC) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)
$(TVL1)      : $(OBJ) $(BINDIR) ; $(CC) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)
$(BINDIR)%.sh: scripts/%.sh $(BINDIR); cp $< $@ 
$(BINDIR)    : ; mkdir -p $@
clean        : ; $(RM) $(NLDCT) $(TVL1) $(OBJ) $(SCRIPT)


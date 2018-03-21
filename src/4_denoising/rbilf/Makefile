CFLAGS   = -std=c99 -O3 -DNDEBUG -w -fopenmp
CXXFLAGS = -O3 -DNDEBUG -w
CPPFLAGS = -Ilib/iio -Ilib/argparse -Ilib/tvl1flow
LDLIBS   = -ljpeg -lpng -ltiff -lm
OMPFLAGS = -fopenmp

BINDIR    = build/bin/
RBILF     = $(BINDIR)rbilf
TVL1      = $(BINDIR)tvl1flow
SCRIPT    = $(BINDIR)rbilf-gt.sh $(BINDIR)rbilf-mp4.sh $(BINDIR)psnr.sh
OBJ_RBILF = lib/iio/iio.o lib/argparse/argparse.o src/main.o
OBJ_TVL1  = lib/iio/iio.o lib/tvl1flow/main.o

all          : $(RBILF) $(TVL1) $(SCRIPT)
$(RNLM)      : $(OBJ_RBILF) $(BINDIR) ; $(CC) $(LDFLAGS) $(OMPFLAGS) -o $@ $(OBJ_RBILF) $(LDLIBS)
$(TVL1)      : $(OBJ_TVL1)  $(BINDIR) ; $(CC) $(LDFLAGS) $(OMPFLAGS) -o $@ $(OBJ_TVL1)  $(LDLIBS)
$(BINDIR)%.sh: scripts/%.sh $(BINDIR); cp $< $@ 
$(BINDIR)    : ; mkdir -p $@
clean        : ; $(RM) $(RBILF) $(TVL1) $(OBJ_TVL1) $(OBJ_RRBILF) $(SCRIPT)


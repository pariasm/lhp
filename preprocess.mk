INPUTS   = $(shell ls i*tif)
LINEAR   = $(addprefix lin_,$(INPUTS))
UNBAND   = $(addprefix  ub_,$(LINEAR))
MASKED   = $(addprefix   m_,$(UNBAND))
FILLED   = $(addprefix   s_,$(MASKED))

all: $(FILLED)

# linearization: contrast change to achieve linear scaling of gray values
lin_i%.tif: i%.tif
	plambda $^ vavg\ log -o $@

# pointwise average of linearized images
avg_lin.tif: $(LINEAR)
	veco -g finite avg $^ -o $@

# pointwise standard deviation of linearized images
std_lin.tif: $(LINEAR)
	veco -g finite std $^ -o $@

# robustified min and max for conversion to and from png
range.txt: avg_lin.tif
	imprintf "%q[1] %q[99]\n" avg_lin.tif > range.txt

# quantized avg_lin (using simplest color balance)
avg_lin.png: avg_lin.tif range.txt
	plambda avg_lin.tif "`cat range.txt` qe" -o avg_lin.png

# remove bands using Yohann's program
ub_avg_lin.png: avg_lin.png
	demo_MIRE $< $@

# extract image of bands
bands.tif: avg_lin.png ub_avg_lin.png range.txt
	plambda ub_avg_lin.png "`cat range.txt` iqe" -o ub_avg_lin.tif
	plambda avg_lin.tif ub_avg_lin.tif - -o bands2d.tif
	veco -c med bands2d.tif | plambda avg_lin.tif - "" -o bands.tif

# remove bands from each frame
ub_lin_i%.tif: lin_i%.tif bands.tif
	plambda $^ - -o $@

# pixelwise standard deviation
std_ub_lin.tif: $(UNBAND)
	veco -g finite std $^ -o $@

# global mask of "bad" pixels
global_mask.png: std_ub_lin.tif
	plambda $< "x x%W5000 > x x%W995000 < and" -o $@

# apply global and local masks
m_ub_%.tif: ub_%.tif global_mask.png
	plambda $^ "x x%W5000 > x x%W995000 < and y and x nan if" -o $@

# inpaint remaining holes
s_m_%.tif: m_%.tif
	simpois -i $< -o $@ 2>/dev/null

# cleanup
clean:
	$(RM) $(LINEAR) avg_lin.tif std_lin.tif ub_avg_lin.png avg_lin.png \
	range.txt bands.tif ub_avg_lin.tif ub_lin_*.tif bands2d.tif \
	m_ub_*.tif s_m_ub_*.tif std_ub_lin.tif global_mask.png

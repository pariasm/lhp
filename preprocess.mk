INPUTS   = $(shell ls i*tif)
LINEAR   = $(addprefix lin_,$(INPUTS))

all: avglin.tif stdlin.tif

# linearization: contrast change to achieve linear scaling of gray values
lin_i%.tif: i%.tif
	plambda $^ vavg\ log -o $@

# pointwise average of linearized images
avglin.tif: $(LINEAR)
	veco -g finite avg $^ -o $@

# pointwise standard deviation of linearized images
stdlin.tif: $(LINEAR)
	veco -g finite std $^ -o $@

# robustified min and max
range.txt: avglin.tif
	imprintf "%q[1] %q[99]\n" $^ > $@

# quantized avglin (using simplest color balance)
avglin.png: avglin.tif range.txt
	plambda avglin.tif "`cat range.txt` qe" -o avglin.png

# find bands using Yohann's program
ub_avglin.png: avglin.png
	demo_MIRE $^ $@

# extract normalized bands
#bands.tif: ub_avglin.png
#	plambda avglin.png ub_avglin.png - | veco -c avg - | plambda av



clean:
	$(RM) $(LINEAR) avglin.tif stdlin.tif ub_avglin.png avglin.png bands.tif

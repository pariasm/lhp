import subprocess
import shutil
import time
from math import floor
import os
import glob

#'gbus',
#'gstennis',
#'gsalesman',
sequences = [
'gbicycle']

base_dir = '~/Remote/avocat/denoising/data/tut/'

#sigma_str  = [10,20,30,40,50]
sigma_str  = [30,40,50]
method = 'oflow4srch'

for iseq in range(0,len(sequences)):
    for sigma in sigma_str:
	seq = sequences[iseq]
	print "Processing " + seq + " with s = {0:d}".format(sigma)

	# time execution
	start_time = time.time()

	# noisy sequence path
	cln_pat = base_dir + seq + '/s00/%03d.png'
	nsy_pat = base_dir + seq + '/s{0:02d}/%03d.tif'.format(sigma)

	# optical flow path
	fof_pat = base_dir + seq + '/s{0:02d}/tvl1_%03d_f.flo'.format(sigma)
	bof_pat = base_dir + seq + '/s{0:02d}/tvl1_%03d_b.flo'.format(sigma)

#	# params 5x4
#	# params 8x2
#	# params 10x1
#	# params 7x4
#	# params 7x2 ~ varying beta
#	px1 = 7
#	pt1 = 2
#	wx1 = 27
#	wt1 = 6
#	np1 = 100
#	bt1 = 0.8
#	rk1 = 39
#
#	px2 = 7
#	pt2 = 2
#	wx2 = 27
#	wt2 = 6
#	bt2 = max(0.8, 2.0 + (1.8 - 2.0)*(sigma - 12.)/(24. - 12.))
#	np2 = max(40, 40 + (sigma - 24.)*(60. - 40.)/(48. - 24.))
#	rk2 = min(39, np2 - 1)

	# params 7x2 ~ beta = 1
	px1 = 7
	pt1 = 2
	wx1 = 27
	wt1 = 6
	np1 = 60
	bt1 = 1.0
	rk1 = 39

	px2 = 7
	pt2 = 2
	wx2 = 27
	wt2 = 6
	np2 = 60
	bt2 = 1.0
	rk2 = min(39, np2 - 1)

	# params 7x2 (rgb) ~ beta = 1
	px1 = 7
	pt1 = 2
	wx1 = 27
	wt1 = 6
	np1 = 100
	bt1 = 1.0
	rk1 = 39

	px2 = 7
	pt2 = 2
	wx2 = 27
	wt2 = 6
	np2 = 60
	bt2 = 1.0
	rk2 = min(39, np2 - 1)

	cmd = ' '.join(['../bin/vnlbayes',
	                '-i', nsy_pat, '-c', cln_pat, '-has-noise',
	                '-f', '1', '-l', '30',
	                '-sigma', str(sigma),
	                '-fof', fof_pat, '-bof', bof_pat,
	                '-px1', str(px1), '-pt1', str(pt1),
	                '-px2', str(px2), '-pt2', str(pt2),
	                '-wx1', str(wx1), '-wt1', str(wt1), 
	                '-wx2', str(wx2), '-wt2', str(wt2),
	                '-r1' , str(rk1), '-np1', str(np1), '-th1', str(th1),
	                '-r2' , str(rk2), '-np2', str(np2), '-th2', str(th2)])

	print cmd

	# launch process
	proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell=True)

	stdout_val = proc.communicate()[0]
	print stdout_val

	# time execution
	end_time = time.time()
	run_time = floor(10*(end_time - start_time))/10

	# append execution time to measures.txt
	fi = open("measures.txt", "ab")
	fi.write('-time = ' + str(run_time) + '\n')
	fi.close()

	# create folder for results
	results_folder = '../../results/' + method + '/table_1p6x3_2p6x3_sptwo/w27x6_' + seq + '_s{0:02d}'.format(sigma)
	if not os.path.exists(results_folder):
		os.makedirs(results_folder)

	# move measures to folder
	shutil.move('measures.txt', results_folder + '/measures')
	shutil.move('measures_frames.txt', results_folder + '/measures_frames')

	# move images to folder
	results_prefixes = ['deno', 'bsic', 'nisy']
	for prefix in results_prefixes:
		files = glob.glob(prefix + '*.png')
		for f in files:
			shutil.move(f, results_folder + '/' + f)

	# remove noisy and difference images
	delete_prefixes = ['diff', 'wei1', 'wei2', 'var2']
	for prefix in delete_prefixes:
		files = glob.glob(prefix + '*.png')
		for f in files:
			os.remove(f)

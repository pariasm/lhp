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

        # ====================
        # GRAYSCALE PARAMETERS
        # ====================

#	# params 5x4
#	px1 = 5
#	pt1 = 4
#	np1 = 200
#	th1 = 2.1
#	px2 = 5
#	pt2 = 4
#	np2 = int(round( 40. + (sigma - 6.)*(60. - 40.)/(48. - 6.) ))
#	th2 = max(0.5, 1.7 + (sigma - 6.)*(1.2 - 1.7)/(24. - 6.))

	# params 7x2
	px1 = 7
	pt1 = 2
	np1 = 150
	th1 = 2.1
	px2 = 7
	pt2 = 2
	np2 = int(round( 40. + (sigma - 6.)*(60. - 40.)/(48. - 6.) ))
	th2 = max(0.5, 2.2 + (sigma - 6.)*(1.2 - 2.2)/(24. - 6.))

#	# params 6x3
#	px1 = 6
#	pt1 = 3
#	np1 = 200
#	th1 = 2.1
#	px2 = 6
#	pt2 = 3
#	np2 = int(round(35. + (sigma - 6.)*(60. - 35.)/(48. - 6.)))
#	th2 = max(0.5,  2.2 + (sigma - 6.)*(1.2 - 2.2)/(24. - 6.))

#	# params 10x1
#	px1 = 10
#	pt1 = 1
#	np1 = 100
#	th1 = 2.9
#	px2 = 10
#	pt2 = 1
#	np2 = int(round( 40. + (sigma - 6.)*(60. - 40.)/(48. - 6.) ))
#	th2 =   max(0.7, 2.2 + (sigma - 6.)*(1.2 - 2.2)/(24. - 6.))

#	# params 6x4
#	# params 7x3
#	# params 9x2
#	# params 12x1
#	# params 8x2
#	px1 = 8
#	pt1 = 2
#	np1 = int(round(max(100.,200. + (sigma - 6.)*(150. - 200.)/(24. - 6.)))) # 200
#	th1 =           min(3.7 ,2.1  + (sigma - 6.)*(2.9  - 2.1 )/(24. - 6.))   # 2.1
#	px2 = 8
#	pt2 = 2
#	np2 = 40
#	th2 = max(0.7, 2.2 + (sigma - 6.)*(1.2 - 2.2)/(24. - 6.))

#	# params 8x3
#	# params 10x2
#	# params 14x1
#	# params 7x4
#	px1 = 7
#	pt1 = 4
#	np1 = 150
#	th1 = 2.1 + (sigma - 6.)*(2.9 - 2.1)/(48. - 6.)
#	px2 = 7
#	pt2 = 4
#	np2 = 60
#	th2 = 1.7 + (sigma - 6.)*(0.7 - 1.7)/(48. - 6.)


	wx1 = 27
	wx2 = 27
	wt1 = 6
	wt2 = 6
	rk1 = min(39, np1 - 1)
	rk2 = min(39, np2 - 1)




        # ================
        # COLOR PARAMETERS
        # ================

        # params 5x4 color
        px1 = 5 
        pt1 = 4 
        np1 = 150 
        th1 = 2.9 
        px2 = 5 
        pt2 = 4 
        np2 = 60
        th2 = max(0.2, 1.7 + (1.2 - 1.7)*(sigma - 10)/(20 - 10))

        # params 7x2 color
        px1 = 7
        pt1 = 2
        np1 = 100
        th1 = 3.7
        px2 = 7
        pt2 = 2
        np2 = 60
        th2 = max(0.2, 1.7 + (1.2 - 1.7)*(sigma - 10)/(20 - 10))

        # params 10x1 color
        px1 = 10
        pt1 = 1
        np1 = 200
        th1 = 2.9
        px2 = 10
        pt2 = 1
        np2 = 80
        th2 = max(0., 0.7 + (0.2 - 0.7)*(sigma - 10)/(20 - 10))



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

/*
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work: Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <sstream>
#include <float.h>

#include "Utilities/Utilities.h"
#include "NlBayes/VideoNLBayes.h"
#include "Utilities/cmd_option.h"

using namespace std;

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author MARC LEBRUN  <marc.lebrun.ik@gmail.com>
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

enum Mode { BSIC_DENO, BSIC_ONLY, DENO_ONLY, NISY_ONLY };

int main(int argc, char **argv)
{
	clo_usage("Video NL-Bayes video denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  clean_path = clo_option("-c"    , ""              , "< clean sequence [to compute PSNR when -has-noise]");
	const string  inbsc_path = clo_option("-b"    , ""              , "< input basic sequence");
	const string  noisy_path = clo_option("-nisy" , "nisy_%03d.png" , "> noisy sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");
	const string   diff_path = clo_option("-diff" , "diff_%03d.png" , "> difference sequence");
	// TODO: these should be determined automatically from the other outputs.
	const string   bias_path = clo_option("-bdeno", "bdeno_%03d.png", "> bias sequence");
	const string bbasic_path = clo_option("-bbsic", "bbsic_%03d.png", "> bias basic sequence");
	const string  bdiff_path = clo_option("-bdiff", "bdiff_%03d.png", "> bias difference sequence");

	      unsigned first_frame = clo_option("-f", 0, "< first frame");
	const unsigned last_frame  = clo_option("-l", 0, "< last frame");
	const unsigned frame_step  = clo_option("-s", 1, "< frame step");

	//! Paths to optical flow
	const string  fflow_path = clo_option("-fof", "", "< input forward  optical flow");
	const string  bflow_path = clo_option("-bof", "", "< input backward optical flow");

	//! General parameters
	const float sigma  = clo_option("-sigma", 0.f, "< add noise of standard deviation sigma");
	const float sigmab = clo_option("-sigma_basic", 0.f, "< standard dev. of remanent noise in the basic estimate");
	      bool has_noise  = (bool) clo_option("-has-noise"   , false, "< input image already has noise");
	const bool do_bias    = (bool) clo_option("-compute-bias", false, "< compute bias outputs");
	const bool verbose    = (bool) clo_option("-verbose"     , true , "< verbose output");
	const bool use_oracle = (bool) clo_option("-oracle"      , false, "< use clean sequence as oracle");
	const unsigned print_prms = (unsigned) clo_option("-print-prms", 0, "< prints parameters for given channels");

	//! Video NLB parameters
	const int time_search1  = clo_option("-wt1", 0  , "< search window temporal radius, step 1");
	const int time_search2  = clo_option("-wt2", 0  , "< search window temporal radius, step 2");
	const int space_search1 = clo_option("-wx1",-1  , "< search window spatial radius, step 1");
	const int space_search2 = clo_option("-wx2",-1  , "< search window spatial radius, step 2");
	const int patch_sizex1  = clo_option("-px1",-1  , "< spatial patch size, step 1");
	const int patch_sizex2  = clo_option("-px2",-1  , "< spatial patch size, step 2");
	const int patch_sizet1  = clo_option("-pt1", 1  , "< temporal patch size, step 1");
	const int patch_sizet2  = clo_option("-pt2", 1  , "< temporal patch size, step 2");
	const int num_patches1  = clo_option("-np1",-1  , "< number of similar patches, step 1");
	const int num_patches2  = clo_option("-np2",-1  , "< number of similar patches, step 2");
	const int rank1         = clo_option("-r1" , 4  , "< rank of covariance matrix, step 1");
	const int rank2         = clo_option("-r2" , 4  , "< rank of covariance matrix, step 2");
#if (defined(CLIPPED_VARIANCE) || defined(PAUL_SIMPLE_VARIANCE) || defined(PAUL_VARIANCE) || defined(VARIANCE_THRESHOLD_1))
	const float thres1      = clo_option("-th1", 0.f, "< variance threshold for cov matrix, step 1");
#else
	const float thres1      = clo_option("-th1", -FLT_MAX, "< variance threshold for cov matrix, step 1");
#endif
	const float thres2      = clo_option("-th2", 0.f, "< variance threshold for cov matrix, step 2");
	const float beta1       = clo_option("-b1" , 1.f, "< noise variance correction factor, step 1");
	const float beta2       = clo_option("-b2" , 1.f, "< noise variance correction factor, step 2");
	const int patch_step1   = clo_option("-sp1",-1  , "< patch skipping step, step 1");
	const int patch_step2   = clo_option("-sp2",-1  , "< patch skipping step, step 2");
	const bool flat_area1 = (bool) clo_option("-flat-area1", false , "< use flat area trick, step 1");
	const bool flat_area2 = (bool) clo_option("-flat-area2", false , "< use flat area trick, step 2");
	const bool no_paste1  = (bool) clo_option("-no-paste1", false , "< disable paste trick, step 1");
	const bool no_paste2  = (bool) clo_option("-no-paste2", false , "< disable paste trick, step 2");
	const int  only_frame = clo_option("-only",-1, "< process only given frame");

	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input sequence.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if (patch_sizex1 == 0 && patch_sizex2 > 0 && inbsc_path == "")
	{
		fprintf(stderr, "%s: if px1 = 0 and px2 > 0, a basic sequence path must be given.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if ((patch_sizex1 < 0 && patch_sizex1 != -1) ||
	    (patch_sizex2 < 0 && patch_sizex2 != -1) ||
	    (patch_sizet1 < 0 || patch_sizet2 <   0) )
	{
		fprintf(stderr, "%s: px1, px2, pt1 and pt2 cannot be negative.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if ((num_patches1 < 0 && num_patches1 != -1) ||
	    (num_patches2 < 0 && num_patches2 != -1) ||
	    (rank1  < 0 || rank2  < 0) ||
#if (defined(CLIPPED_VARIANCE) || defined(PAUL_SIMPLE_VARIANCE))
	    // thres1 can be negative if negative eigenvalues are used in the
	    // second step. Otherwise, it has to be positive.
	    (thres1 < 0) ||
#endif
	    (thres2 < 0) )
	{
		fprintf(stderr, "%s: np1, np2, r1, r2, th1 and th2 cannot be negative.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if ((space_search1 < 0 && space_search1 != -1) ||
	    (space_search2 < 0 && space_search2 != -1) ||
	    ( time_search1 < 0 ||  time_search2 <   0) )
	{
		fprintf(stderr, "%s: wx1, wx2, wt1 and wt2 cannot be negative.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if (patch_sizex1 > 0 && inbsc_path != "")
		fprintf(stderr, "\x1b[33;1mWarning:\x1b[0m basic sequence path ignored since px1 > 0.\n");

	if ((fflow_path != "" && bflow_path == "") || 
	    (fflow_path == "" && bflow_path != ""))
	{
		fprintf(stderr, "Only one oflow path provided.\nTry `%s --help' for more information.\n",
				argv[0]);
		return EXIT_FAILURE;
	}

	if (use_oracle && has_noise && clean_path == "")
	{
		fprintf(stderr, "The clean video must be provided in oracle mode.\nTry `%s --help' for more information.\n",
				argv[0]);
		return EXIT_FAILURE;
	}

	//! Determine mode
	Mode mode;
	if ((patch_sizex1 != 0) && (patch_sizex2 != 0)) mode = BSIC_DENO;
	if ((patch_sizex1 != 0) && (patch_sizex2 == 0)) mode = BSIC_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 != 0)) mode = DENO_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 == 0)) mode = NISY_ONLY;

	if (clean_path != "") has_noise = true;

	bool use_oflow = (fflow_path != "");

	//! Only print parameters
	if (print_prms)
	{
		VideoSize tmp;
		tmp.channels = print_prms;

		//! Compute denoising default parameters
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, sigma, tmp, flat_area1, verbose, time_search1, time_search1, patch_sizet1);
		VideoNLB::initializeNlbParameters(prms2, 2, sigma, tmp, flat_area2, verbose, time_search2, time_search2, patch_sizet2);

		//! Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, tmp, (unsigned)patch_sizex1);;
		if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, tmp, (unsigned)patch_sizex2);;
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);

		int pdim1 = prms1.sizePatch * prms1.sizePatch * prms1.sizePatchTime;
#ifdef CHANNELS_DECOUPLED
		int pdim2 = prms2.sizePatch * prms2.sizePatch * prms2.sizePatchTime;
#else
		int pdim2 = prms2.sizePatch * prms2.sizePatch * prms2.sizePatchTime * tmp.channels;
#endif

//		prms1.rank = (thres1 == 0.f) ? rank1 : pdim1 - 1;
//		prms2.rank = (thres2 == 0.f) ? rank2 : pdim2 - 1;
		prms1.rank = rank1;
		prms2.rank = rank2;

		prms1.variThres = thres1;
		prms2.variThres = thres2;

		prms1.beta = beta1;
		prms2.beta = beta2;

		prms1.onlyFrame = only_frame - first_frame;
		prms2.onlyFrame = only_frame - first_frame;

		if (no_paste1) prms1.doPasteBoost = false;
		if (no_paste2) prms2.doPasteBoost = false;

		if (patch_step1 >= 0) prms1.offSet = patch_step1;
		if (patch_step2 >= 0) prms2.offSet = patch_step2;

		prms2.sigma_basic = sigmab;
 
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		return EXIT_SUCCESS;
	}


	//! Declarations
	Video<float> original, noisy, basic, final, diff;
	Video<float> fflow, bflow;

	//! Load input videos
	                       original.loadVideo(input_path, first_frame, last_frame, frame_step);
	if (mode == DENO_ONLY) basic   .loadVideo(inbsc_path, first_frame, last_frame, frame_step);
	if (use_oflow)         fflow   .loadVideo(fflow_path, first_frame, last_frame, frame_step);
	if (use_oflow)         bflow   .loadVideo(bflow_path, first_frame, last_frame, frame_step);

	//! Add noise
	if (has_noise)
	{
		if (verbose) printf("Input video has noise of sigma = %f\n", sigma);
		noisy = original;

		if (clean_path != "")
			original.loadVideo(clean_path, first_frame, last_frame, frame_step);
		else 
			if (verbose) printf("No clean video provided.\n");
	}
	else if (sigma)
		VideoUtils::addNoise(original, noisy, sigma, verbose);

	//! Save noisy video
	if (mode == NISY_ONLY)
	{
		if (!has_noise && sigma)
		{
			float noisy_psnr = -1, noisy_rmse = -1;
			VideoUtils::computePSNR(original, noisy, noisy_psnr, noisy_rmse);
			writingMeasures("measures.txt", sigma, noisy_psnr, noisy_rmse, 0, true, "_noisy");

			if (verbose) printf("Saving noisy video\n");
			noisy.saveVideo(noisy_path, first_frame, frame_step);

			return EXIT_SUCCESS;
		}
		else
		{
			fprintf(stderr, "Noisy video not saved since is equal to original input.\n");
			return EXIT_FAILURE;
		}
	}

	//! Denoising
	if (verbose) printf("Running Video NL-Bayes on the noisy video\n");
	if (verbose && use_oracle) printf("Using provided oracle\n");

	//! Compute denoising default parameters
	VideoNLB::nlbParams prms1, prms2;
	VideoNLB::initializeNlbParameters(prms1, 1, sigma, noisy.sz, flat_area1, verbose, time_search1, time_search1, patch_sizet1);
	VideoNLB::initializeNlbParameters(prms2, 2, sigma, noisy.sz, flat_area2, verbose, time_search2, time_search2, patch_sizet2);

	//! Override with command line parameters
	if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
	if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
	if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, noisy.sz, (unsigned)patch_sizex1);;
	if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, noisy.sz, (unsigned)patch_sizex2);;
	if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
	if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);

	int pdim1 = prms1.sizePatch * prms1.sizePatch * prms1.sizePatchTime;
#ifdef CHANNELS_DECOUPLED
	int pdim2 = prms2.sizePatch * prms2.sizePatch * prms2.sizePatchTime;
#else
	int pdim2 = prms2.sizePatch * prms2.sizePatch * prms2.sizePatchTime * noisy.sz.channels;
#endif

//	prms1.rank = (thres1 == 0.f) ? rank1 : pdim1 - 1;
//	prms2.rank = (thres2 == 0.f) ? rank2 : pdim2 - 1;
	prms1.rank = rank1;
	prms2.rank = rank2;

	prms1.variThres = thres1;
	prms2.variThres = thres2;

	prms1.beta = beta1;
	prms2.beta = beta2;

	if (no_paste1) prms1.doPasteBoost = false;
	if (no_paste2) prms2.doPasteBoost = false;

	if (patch_step1 >= 0) prms1.offSet = patch_step1;
	if (patch_step2 >= 0) prms2.offSet = patch_step2;

	prms1.onlyFrame = only_frame - first_frame;
	prms2.onlyFrame = only_frame - first_frame;

	prms2.sigma_basic = sigmab;

	//! Percentage or processed groups of patches over total number of pixels
	std::vector<float> groupsRatio;

	//! Oracle mode
	Video<float> oracle = (use_oracle) ? original : Video<float>();

	//! Run denoising algorithm
	if (use_oflow)
		groupsRatio = VideoNLB::runNlBayes(noisy, fflow, bflow, basic, final, prms1, prms2, oracle);
	else
		groupsRatio = VideoNLB::runNlBayes(noisy, basic, final, prms1, prms2, oracle);

	if (only_frame >= 0)
	{
		Video<float> tmp(noisy.sz.width, noisy.sz.height, 1, noisy.sz.channels);
		if (prms2.sizePatch)
		{
			VideoUtils::crop(final, tmp, only_frame - first_frame);
			final = tmp;
		}

		VideoUtils::crop(basic, tmp, only_frame - first_frame);
		basic = tmp;

		VideoUtils::crop(noisy, tmp, only_frame - first_frame);
		noisy = tmp;

		VideoUtils::crop(original, tmp, only_frame - first_frame);
		original = tmp;

		first_frame = only_frame;
	}

	//! Compute PSNR and RMSE
	float final_psnr = -1, final_rmse = -1, basic_psnr = -1, basic_rmse = -1;
	if (prms2.sizePatch) VideoUtils::computePSNR(original, final, final_psnr, final_rmse);
	                     VideoUtils::computePSNR(original, basic, basic_psnr, basic_rmse);

	if (verbose)
	{
	    printf("basic PSNR =\t%f\tRMSE =\t%f\n", basic_psnr, basic_rmse);
	    printf("final PSNR =\t%f\tRMSE =\t%f\n", final_psnr, final_rmse);
	}

	//! Write measures
	writingMeasures("measures.txt", sigma, basic_psnr, basic_rmse, groupsRatio[0], true  , "_basic");
	writingMeasures("measures.txt", sigma, final_psnr, final_rmse, groupsRatio[1], false , "_final");

	//! Per-frame measures
	std::vector<float> final_frames_psnr, final_frames_rmse, basic_frames_psnr, basic_frames_rmse;
	if (prms2.sizePatch) VideoUtils::computePSNR(original, final, final_frames_psnr, final_frames_rmse);
	                     VideoUtils::computePSNR(original, basic, basic_frames_psnr, basic_frames_rmse);
	writingMeasures("measures_frames.txt", sigma, basic_frames_psnr, basic_frames_rmse, groupsRatio[0], true  , "_basic");
	writingMeasures("measures_frames.txt", sigma, final_frames_psnr, final_frames_rmse, groupsRatio[1], false , "_final");

	//! Compute Difference
	if (prms2.sizePatch) VideoUtils::computeDiff(original, final, diff, sigma);

	//! Save output sequences
	if (verbose) printf("Saving output sequences\n");

	if (prms2.sizePatch) final.saveVideo(final_path, first_frame, frame_step);
	if (prms2.sizePatch) diff .saveVideo( diff_path, first_frame, frame_step);
	                     basic.saveVideo(basic_path, first_frame, frame_step);

#if 0
	//! Computing bias sequence
	Video<float> bias, bias_basic, bias_diff;
	if (do_bias) {
		if (verbose) cout << "Applying NL-Bayes to the original image :" << endl;

		runNlBayes(im, imBasicBias, imBias, imSize, flat_area1, flat_area2, sigma, verbose);

		if (verbose) cout << endl;

		float psnrBias, psnrBiasBasic, rmseBias, rmseBiasBasic;
		computePsnr(im, imBasicBias, psnrBiasBasic, rmseBiasBasic, "imBiasBasic", verbose);
		computePsnr(im, imBias     , psnrBias     , rmseBias     , "imBiasFinal", verbose);

		//! writing measures
		writingMeasures("measures.txt", sigma, psnrBiasBasic, rmseBiasBasic, false, "_bias_basic");
		writingMeasures("measures.txt", sigma, psnrBias     , rmseBias     , false, "_bias      ");

		
		computeDiff(im, imBias, imDiffBias, sigma, 0.f, 255.f, verbose);

		saveImage(argv[7], imBias     , imSize, 0.f, 255.f);
		saveImage(argv[8], imBasicBias, imSize, 0.f, 255.f);
		saveImage(argv[9], imDiffBias , imSize, 0.f, 255.f);
	}


#endif

	if (verbose) printf("Done\n");
	return EXIT_SUCCESS;
}

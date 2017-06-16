/*
 * Modified work Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file VideoNLBayes.cpp
 * @brief Video NL-Bayes denoising functions
 *
 * @author Pablo Arias <pariasm@gmail.com>
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <float.h>

#include "VideoNLBayes.h"
#include "LibMatrix.h"
#include "../Utilities/Utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Choose implementation for low-rank Bayes estimate in both steps. If both
 * options are commented, the spectral decomposition of the covariance matrix
 * is computed using LAPACK's low-rank eigs */
//#define USE_SVD_LAPACK

// search region centered on the reference patch; image extended by symmetrizing
//#define CENTRED_SEARCH

// dumps aggregation weights (also useful to dump other scalar quantities)
//#define DEBUG_SHOW_WEIGHT

// returns an image with only color coded reference patches
//#define DEBUG_SHOW_REFERENCE_PATCHES

// dumps all the non-flat patches, noisy and after Bayesian estimation
//#define DEBUG_SHOW_GROUPS 


// colors
#define ANSI_BLK  "\x1b[30m"
#define ANSI_RED  "\x1b[31m"
#define ANSI_GRN  "\x1b[32m"
#define ANSI_YLW  "\x1b[33m"
#define ANSI_BLU  "\x1b[34m"
#define ANSI_MAG  "\x1b[35m"
#define ANSI_CYN  "\x1b[36m"
#define ANSI_WHT  "\x1b[37m"
#define ANSI_BBLK "\x1b[30;01m"
#define ANSI_BRED "\x1b[31;01m"
#define ANSI_BGRN "\x1b[32;01m"
#define ANSI_BYLW "\x1b[33;01m"
#define ANSI_BBLU "\x1b[34;01m"
#define ANSI_BMAG "\x1b[35;01m"
#define ANSI_BCYN "\x1b[36;01m"
#define ANSI_BWHT "\x1b[37;01m"
#define ANSI_RST  "\x1b[0m"

namespace VideoNLB
{

/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param o_params   : will contain the nlbParams for the first step of the algorithm;
 * @param p_step     : select first or second step;
 * @param p_sigma    : standard deviation of the noise;
 * @param p_size     : size of the video;
 * @param p_flatArea : if true, use the homogeneous area trick for the first step;
 * @param p_verbose  : if true, print some informations.
 * @param p_timeSearchRagneFwd : temporal search range forwards.
 * @param p_timeSearchRagneBwd : temporal search range backwards.
 *
 * @return none.
 **/
void initializeNlbParameters(
	nlbParams &o_params
,	const unsigned p_step
,	const float p_sigma
,	const VideoSize &p_size
,	const bool p_flatArea
,	const bool p_verbose
,	const unsigned timeSearchRangeFwd
,	const unsigned timeSearchRangeBwd
,	const unsigned sizePatchTime
,	const unsigned rank
){
	const bool s1 = (p_step == 1);

	//! Standard deviation of the noise
	o_params.sigma = p_sigma;

	//! Patch size, spatial dimension
	if (p_size.channels == 1)
	{
		if(s1) o_params.sizePatch = (p_sigma < 30.f ? 5 : 7);
		else   o_params.sizePatch = 5;
	}
	else
	{
		if(s1) o_params.sizePatch = (p_sigma < 20.f ? 3 : (p_sigma < 50.f ? 5 : 7));
		else   o_params.sizePatch = (p_sigma < 50.f ? 3 : (p_sigma < 70.f ? 5 : 7));
	}

	//! Patch size, temporal dimension
	o_params.sizePatchTime = sizePatchTime;

	//! Number of similar patches
	if (p_size.channels == 1)
	{
		if(s1) o_params.nSimilarPatches = (p_sigma < 10.f ?  35 : 
		                                  (p_sigma < 30.f ?  45 : 
		                                  (p_sigma < 80.f ?  90 : 
		                                                    100))); // ASK MARC differs from manuscript
		else   o_params.nSimilarPatches = (p_sigma < 20.f ?  15 : 
		                                  (p_sigma < 40.f ?  25 :
		                                  (p_sigma < 80.f ?  30 :
		                                                     45))); // ASK MARC differs from manuscript
	}
	else
		o_params.nSimilarPatches = o_params.sizePatch * o_params.sizePatch * 3;

	//! Offset: step between two similar patches
	o_params.offSet     = std::max((unsigned)1, o_params.sizePatch     / 2);
//	o_params.offSetTime = std::max((unsigned)1, o_params.sizePatchTime / 2);
	o_params.offSetTime = 1;

	//! Use the homogeneous area detection trick
	o_params.useHomogeneousArea = p_flatArea;

	//! Size of the search window around the reference patch (must be odd)
	o_params.sizeSearchWindow = o_params.nSimilarPatches / 2; // ASK MARC and 7*sizePatch1 in IPOL manuscript
	if (o_params.sizeSearchWindow % 2 == 0) o_params.sizeSearchWindow++;

	//! Search window, temporal search radii
	o_params.sizeSearchTimeRangeFwd = timeSearchRangeFwd;
	o_params.sizeSearchTimeRangeBwd = timeSearchRangeBwd;
	o_params.nSimilarPatches *= timeSearchRangeFwd + timeSearchRangeBwd + 1;

	//! Size of boundaries used during the sub division
	o_params.boundary = 2*(o_params.sizeSearchWindow/2) + (o_params.sizePatch - 1);

	//! Parameter used to determine if an area is homogeneous
#ifdef FAT_ORIGINAL
	o_params.gamma = 1.05f;
#else
	o_params.gamma = 0.95f;
#endif

	//! Variance threshold on principal directions
	o_params.variThres = .0f;

	//! Parameter used to estimate the covariance matrix
	if (p_size.channels == 1)
	{
		o_params.beta = 1.f;
//		if(s1) o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 70.f ? 1.f : 0.9f));
//		else   o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 35.f ? 1.f : 0.9f));
	}
	else
	{
		o_params.beta = 1.f;
//		if(s1) o_params.beta = 1.f;
//		else   o_params.beta = (p_sigma < 50.f ? 1.2f : 1.f);
	}

	// maximum rank of covariance matrix
	o_params.rank = rank;

	//! Parameter used to determine similar patches
	//  Differs from manuscript (tau = 4) because (1) this threshold is to be used 
	//  with the squared distance and (2) the distance is not normalized
#ifndef SPTWO
	if(s1) o_params.tau = 0; // not used
	else   o_params.tau = 16.f * o_params.sizePatch * o_params.sizePatch * p_size.channels;
#else
	o_params.tau = 0;
#endif

	//! Print information?
	o_params.verbose = p_verbose;

	//! Is first step?
	o_params.isFirstStep = (p_step == 1);

	//! Boost the paste trick
	o_params.doPasteBoost = true;

	//! Process all frames
	o_params.onlyFrame = -1;
}

//	depend on sigma:
//		  sizePatch
//		  nSimilarPatches
//		  beta
//
// depend on sizePatch:
//		  nSimilarPatches (if channels == 3)
//		  sizeSearchWindow
//		  offSet
//		  tau
//
// depend on nSimilarPatches
//		  sizeSearchWindow
//
// depend on sizeSearchWindow
//		  boundary
//		  nSimilarPatches (cannot be more than total number of searchable patches)
//
// depend on sizeSearchTimeRangeFwd/Bwd
//		  nSimilarPatches

/**
 * @brief Sets size of spatial search window. It sets the border width accordingly,
 * and also ensures that the number of similar patches is not larger that the 
 * total number of available patches.
 *
 * @param prms             : nlbParams for first or second step of the algorithm;
 * @param sizeSearchWindow : size of search window;
 *
 * @return none.
 **/
void setSizeSearchWindow(nlbParams& prms, unsigned sizeSearchWindow)
{
	prms.sizeSearchWindow = sizeSearchWindow;
	prms.boundary = 2*(sizeSearchWindow/2) + (prms.sizePatch - 1);
	prms.nSimilarPatches = std::min(prms.nSimilarPatches, sizeSearchWindow *
	                                                      sizeSearchWindow *
	                                                     (prms.sizeSearchTimeRangeFwd +
	                                                      prms.sizeSearchTimeRangeBwd + 1));
}

/**
 * @brief Sets size of the patch. It sets the pixel offset as half the patch
 * size (this is BM3D speed-up).
 *
 * @param prms      : nlbParams for first or second step of the algorithm;
 * @param sizePatch : size of the patch;
 *
 * @return none.
 **/
void setSizePatch(nlbParams& prms, const VideoSize &size, unsigned sizePatch)
{
	prms.sizePatch = sizePatch;
	prms.boundary = 2*(prms.sizeSearchWindow/2) + (prms.sizePatch - 1);
	prms.offSet = sizePatch/2;

	//! Update number of similar patches, only if it is less than recommended value
	if (size.channels == 3)
//		prms.nSimilarPatches = std::max(prms.sizePatch * prms.sizePatch * 3, prms.nSimilarPatches);
		prms.nSimilarPatches = prms.sizePatch * prms.sizePatch * 3 * 
		                      (prms.sizeSearchTimeRangeFwd + 
		                       prms.sizeSearchTimeRangeBwd + 1);
}

/**
 * @brief Sets number of similar patches, ensuring that the number of similar
 * patches is not larger that the total number of available patches.
 *
 * @param prms            : nlbParams for first or second step of the algorithm;
 * @param nSimilarPatches : number of similar patches;
 *
 * @return none.
 **/
void setNSimilarPatches(nlbParams& prms, unsigned nSimilarPatches)
{
	prms.nSimilarPatches = nSimilarPatches;
	prms.nSimilarPatches = std::min(nSimilarPatches, prms.sizeSearchWindow *
	                                                 prms.sizeSearchWindow *
	                                                (prms.sizeSearchTimeRangeFwd +
	                                                 prms.sizeSearchTimeRangeBwd + 1));
}

/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_params : nlbParams for first or second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	const nlbParams &i_prms
){
//	printf("\nVideo NLBayes parameters\n");
//	printf("------------------------\n\n");
//	printf("Noise sigma = %g\n", i_prms1.sigma);

	printf("\x1b[37;01m" "Parameters for step %d:" ANSI_RST "\n" , i_prms.isFirstStep ? 1 : 2);
	printf("\tPatch search:\n");
	printf("\t\tPatch size                  = %d\n"       , i_prms.sizePatch);
	printf("\t\tPatch size temporal         = %d\n"       , i_prms.sizePatchTime);
	printf("\t\tNumber of patches           = %d\n"       , i_prms.nSimilarPatches);
	if (!i_prms.isFirstStep) printf("\t\tDistance threshold (tau)    = %g\n"       , i_prms.tau);
	else                     printf("\t\tDistance threshold (tau)    = N/A\n"      );
	printf("\t\tSpatial search window       = %dx%d\n"    , i_prms.sizeSearchWindow, i_prms.sizeSearchWindow);
	printf("\t\tTemporal search range       = [-%d,%d]\n" , i_prms.sizeSearchTimeRangeBwd, i_prms.sizeSearchTimeRangeBwd);
	printf("\t\tSpatial border added        = %d\n"       , i_prms.boundary);
	printf("\tGroup filtering:\n");
	printf("\t\tBeta                        = %g\n"       , i_prms.beta);
	printf("\t\tRank                        = %d\n"       , i_prms.rank);
	printf("\t\tVariance th.                = %g\n"       , i_prms.variThres);
	if (!i_prms.isFirstStep)
		printf("\t\tSigma basic                 = %g\n"       , i_prms.sigma_basic);
	if (i_prms.useHomogeneousArea)
		printf("\t\tFlat area trick with gamma  = %g\n"       , i_prms.gamma);
	else
		printf("\t\tFlat area trick             = inactive\n");
	printf("\tSpeed-ups:\n");
	printf("\t\tOffset                      = %d\n"       , i_prms.offSet);
	printf("\t\tOffsetTime                  = %d\n"       , i_prms.offSetTime);
	printf("\t\tPasteBoost                  = %s\n\n"     , i_prms.doPasteBoost ? "active" : "inactive");
	if (i_prms.onlyFrame >= 0)
		printf("\tProcessing only frame %d\n\n", i_prms.onlyFrame);
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_noisy: contains the noisy image;
 * @param o_basic: will contain the basic estimate image after the first step;
 * @param o_final: will contain the final denoised image after the second step;
 * @param p_useArea1 : if true, use the flat area trick for the first step;
 * @param p_useArea2 : if true, use the flat area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> & i_noisy
,	Video<float> &o_basic
,	Video<float> &o_final
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,	const bool p_verbose
,  Video<float> &i_imClean
){
	//! Video size
	VideoSize size = i_noisy.sz;

	//! Parameters Initialization
	nlbParams p_prms1, p_prms2;
	initializeNlbParameters(p_prms1, 1, p_sigma, size, p_useArea1, p_verbose);
	initializeNlbParameters(p_prms2, 2, p_sigma, size, p_useArea2, p_verbose);

	//! NL-Bayes
	return runNlBayes(i_noisy, o_basic, o_final, p_prms1, p_prms2, i_imClean);
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_prms1  : parameters for first  step;
 * @param p_prms2  : parameters for second step;
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> &i_imNoisy
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_prms1
,	const nlbParams p_prms2
,  Video<float> &i_imClean
){
	return runNlBayes(i_imNoisy, Video<float>(), Video<float>(), 
			o_imBasic, o_imFinal, p_prms1, p_prms2, i_imClean);
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param i_fflow  : forward optical flow;
 * @param i_bflow  : backward optical flow;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_prms1  : parameters for first  step;
 * @param p_prms2  : parameters for second step;
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> & i_imNoisy
,  Video<float> const& i_fflow
,  Video<float> const& i_bflow
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_prms1
,	const nlbParams p_prms2
,  Video<float> &i_imClean
){
	//! Only 1, 3 or 4-channels images can be processed.
	const unsigned chnls = i_imNoisy.sz.channels;
	if (! (chnls == 1 || chnls == 3 || chnls == 4))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong number of "
				"channels. Must be 1, 3 or 4.");

	//! Check if optical flow is valid
	const bool use_oflow = (i_fflow.sz.width > 0);
	if (use_oflow && 
		(i_fflow.sz.channels != 2 || 
		 i_fflow.sz.width  != i_imNoisy.sz.width  ||
		 i_fflow.sz.height != i_imNoisy.sz.height ||
		 i_fflow.sz.frames != i_imNoisy.sz.frames))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong size of fwd optical flow.");

	if (use_oflow && 
		(i_bflow.sz.channels != 2 || 
		 i_bflow.sz.width  != i_imNoisy.sz.width  ||
		 i_bflow.sz.height != i_imNoisy.sz.height ||
		 i_bflow.sz.frames != i_imNoisy.sz.frames))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong size of bwd optical flow.");

	//! Print compiler options
	if (p_prms1.verbose)
	{
#ifdef NOISY_COVARIANCE
		printf(ANSI_BCYN "NOISY_COVARIANCE > Computing 2nd step cov. matrix from noisy patches.\n" ANSI_RST);
#endif
#ifdef LI_ZHANG_DAI
		printf(ANSI_BCYN "LI_ZHANG_DAI > Using Li-Zhang-Dai's signal estimate for empirical Wiener.\n" ANSI_RST);
#endif
#ifdef BARICENTER_BASIC
		printf(ANSI_BCYN "BARICENTER_BASIC > Centering noisy patches with basic baricenter.\n" ANSI_RST);
#endif
#ifdef CLIPPED_VARIANCE
		printf(ANSI_BCYN "CLIPPED_VARIANCE > Clipping negative variances\n" ANSI_RST);
#endif
#ifdef PAUL_VARIANCE
		printf(ANSI_BCYN "PAUL_VARIANCE > Full inverse of Paul expected eigenvalue\n" ANSI_RST);
#endif
#ifdef FRAMES_DECOUPLED
		printf(ANSI_BCYN "PAUL_SIMPLE_VARIANCE > Simplified inverse of Paul expected eigenvalue\n" ANSI_RST);
#endif
#ifdef CHANNELS_DECOUPLED
		printf(ANSI_BCYN "CHANNELS_DECOUPLED > Independent Gaussian models for each channel (2nd step)\n" ANSI_RST);
#endif
#ifdef SPTWO
		printf(ANSI_BCYN "SPTWO > 3D motion comp. patches for search, but 2D for MAP denoising\n" ANSI_RST);
#endif
	}


	//! Set-up multithreading
#ifndef DEBUG_SHOW_GROUPS
	unsigned nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	if (p_prms1.verbose) printf(ANSI_CYN "OpenMP is using %d threads\n" ANSI_RST, nThreads);
#endif
	const unsigned nParts = 2 * nThreads;
#else 
	//! Force one thread
	/* In doing so the patch coordinates in the dumped patch files
	/* are with respect to a unique image coordinate. */
	unsigned nThreads = 1;
	if (p_prms1.verbose) printf(ANSI_CYN "FORCING 1 THREAD\n" ANSI_RST);
	const unsigned nParts = 1;
#endif

	//! Video size
	VideoSize imSize = i_imNoisy.sz;

	//! Determine steps
	int step = p_prms1.sizePatch ?
	          (p_prms2.sizePatch ? 0 : 1) :
	          (p_prms2.sizePatch ? 2 :-1) ;

	if (step == -1)
		throw std::runtime_error("VideoNLB::runNlBayes: Both patch sizes are zero.");

	//! Initialization
	if (step != 2) o_imBasic.resize(imSize);
	if (step != 1) o_imFinal.resize(imSize);

	if (step == 2 && imSize != o_imBasic.sz) 
		throw std::runtime_error("VideoNLB::runNlBayes: sizes of noisy and "
				"basic videos don't match");

	//! Print parameters
	if (p_prms1.verbose) if (step != 2) printNlbParameters(p_prms1);
	if (p_prms2.verbose) if (step != 1) printNlbParameters(p_prms2);

	//! Percentage of processed groups over total number of pixels
	std::vector<float> groupsRatio(2,0.f);

	//! Step 1
	if (step != 2)
	{
		if (p_prms1.verbose)
		{
			printf("1st Step\n");
			for (int p = 0; p < nParts; ++p) printf("\n");
		}

		//! RGB to YUV
		VideoUtils::transformColorSpace(i_imNoisy, true);
		VideoUtils::transformColorSpace(i_imClean, true);

		//! Divide the noisy image into sub-images in order to easier parallelize the process
		std::vector<Video<float> > imNoisySub(nParts);
		std::vector<VideoUtils::CropPosition > imCrops(nParts);
		VideoUtils::subDivideTight(i_imNoisy, imNoisySub, imCrops, p_prms1.boundary, nParts);

		//! Divide oracle
		std::vector<Video<float> > imCleanSub(nParts);
		if (i_imClean.sz.whcf)
			VideoUtils::subDivideTight(i_imClean, imCleanSub, imCrops, p_prms1.boundary, nParts);

		//! Divide optical flow 
		std::vector<Video<float> > fflowSub(nParts), bflowSub(nParts);
		if (use_oflow)
		{
			std::vector<VideoUtils::CropPosition > oflowCrops(nParts);
			VideoUtils::subDivideTight(i_fflow, fflowSub, oflowCrops, p_prms1.boundary, nParts);
			VideoUtils::subDivideTight(i_bflow, bflowSub, oflowCrops, p_prms1.boundary, nParts);
		}

		//! Process all sub-images
		std::vector<Video<float> > imBasicSub(nParts);
		std::vector<Video<float> > imFinalSub(nParts);
		std::vector<unsigned> groupsProcessedSub(nParts);
#ifdef _OPENMP
		// we make a copy of prms structure because, since it is constant,
		// it causes a compilation error with OpenMP (only on IPOL server)
		nlbParams prms1(p_prms1);
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (prms1)
#endif
		for (int n = 0; n < (int)nParts; n++)
			groupsProcessedSub[n] = 
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], p_prms1, imCrops[n],
				               imCleanSub[n]);

		//! Get the basic estimate
		VideoUtils::subBuildTight(imBasicSub, o_imBasic, p_prms1.boundary);

		//! YUV to RGB
		VideoUtils::transformColorSpace(o_imBasic, false);
		VideoUtils::transformColorSpace(i_imNoisy, false);
		VideoUtils::transformColorSpace(i_imClean, false);

		for (int n = 0; n < (int)nParts; n++)
			groupsRatio[0] += 100.f * (float)groupsProcessedSub[n]/(float)imSize.whf;

#ifdef DEBUG_SHOW_WEIGHT
		{
			std::vector<Video<float> > subWeights(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/weight_step1_%d.%d.%d_%%03d.png", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subWeights[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> weight(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subWeights, weight, p_prms1.boundary);

			// Write to disk
			weight.saveVideo("wei1_%03d.png", 1);
		}
#endif

	}

	//! Step 2
	if (step != 1)
	{
		if (p_prms2.verbose)
		{
			if (step == 2) for (int p = 0; p <= nParts; ++p) printf("\n");
			printf("\x1b[%dF2nd Step\n",nParts+1);
			for (int p = 0; p < nParts; ++p) printf("\x1b[2K\n");
		}

		//! RGB to YUV
		VideoUtils::transformColorSpace(o_imBasic, true);
		VideoUtils::transformColorSpace(i_imNoisy, true);
		VideoUtils::transformColorSpace(i_imClean, true);

		//! Divide the noisy and basic images into sub-images in order to easier parallelize the process
		std::vector<Video<float> > imNoisySub(nParts);
		std::vector<Video<float> > imBasicSub(nParts);
		std::vector<Video<float> > imCleanSub(nParts);
		std::vector<VideoUtils::CropPosition > imCrops(nParts);

		VideoUtils::subDivideTight(i_imNoisy, imNoisySub, imCrops, p_prms2.boundary, nParts);
		VideoUtils::subDivideTight(o_imBasic, imBasicSub, imCrops, p_prms2.boundary, nParts);

		//! Divide oracle
		if (i_imClean.sz.whcf)
			VideoUtils::subDivideTight(i_imClean, imCleanSub, imCrops, p_prms2.boundary, nParts);

		//! Divide optical flow 
		std::vector<Video<float> > fflowSub(nParts), bflowSub(nParts);
		if (use_oflow)
		{
			std::vector<VideoUtils::CropPosition > oflowCrops(nParts);
			VideoUtils::subDivideTight(i_fflow, fflowSub, oflowCrops, p_prms2.boundary, nParts);
			VideoUtils::subDivideTight(i_bflow, bflowSub, oflowCrops, p_prms2.boundary, nParts);
		}

		//! Process all sub-images
		std::vector<Video<float> > imFinalSub(nParts);
		std::vector<float> groupsProcessedSub(nParts);
#ifdef _OPENMP
		// we make a copy of prms structure because, since it is constant,
		// it causes a compilation error with OpenMP (only on IPOL server)
		nlbParams prms2(p_prms2);
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (prms2)
#endif
		for (int n = 0; n < (int) nParts; n++)
			groupsProcessedSub[n] = 
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], p_prms2, imCrops[n],
				               imCleanSub[n]);

		//! Get the final result
		VideoUtils::subBuildTight(imFinalSub, o_imFinal, p_prms2.boundary);

		for (int n = 0; n < (int)nParts; n++)
			groupsRatio[1] += 100.f * (float)groupsProcessedSub[n]/(float)imSize.whf;

		//! YUV to RGB
		VideoUtils::transformColorSpace(i_imNoisy, false);
		VideoUtils::transformColorSpace(o_imBasic, false);
		VideoUtils::transformColorSpace(o_imFinal, false);
		VideoUtils::transformColorSpace(i_imClean, false);

#ifdef DEBUG_SHOW_WEIGHT
		{
			std::vector<Video<float> > subWeights(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/weight_step2_%d.%d.%d_%%03d.png", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subWeights[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> weight(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subWeights, weight, p_prms2.boundary);

			// Write to disk
			weight.saveVideo("wei2_%03d.png", 1);
		}
		{
			std::vector<Video<float> > subVars(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/var_step2_%d.%d.%d_%%03d.png", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subVars[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> variance(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subVars, variance, p_prms2.boundary);


			// Write to disk
			variance.saveVideo("var2_%03d.png", 1);
			variance.saveVideoAscii("var2", 1);
		}
#endif
		if (p_prms2.verbose) printf("\n");
	}

	return groupsRatio;
}

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: see nlbParams.
 *
 * @return Number of processed groups of similar patches.
 **/
unsigned processNlBayes(
	Video<float> const& i_imNoisy
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	Video<float> &io_imBasic
,	Video<float> &o_imFinal
,	nlbParams const& p_params
,	VideoUtils::CropPosition p_crop
,  Video<float> const &i_imClean
){
	using std::vector;

	//! Parameters initialization
	const bool step1 = p_params.isFirstStep;
	const unsigned sWx = p_params.sizeSearchWindow;
	const unsigned sWt = p_params.sizeSearchTimeRangeFwd +
	                     p_params.sizeSearchTimeRangeBwd + 1;// VIDEO
	const unsigned sPx = p_params.sizePatch;
	const unsigned sPt = p_params.sizePatchTime;
	const VideoSize sz = i_imNoisy.sz;

	unsigned nInverseFailed = 0;
	const float threshold = p_params.sigma * p_params.sigma * p_params.gamma ;
//	                       (p_params.isFirstStep ? i_imNoisy.sz.channels : 1.f);

	//! Weight sum per pixel
	Video<float> weight(sz.width, sz.height, sz.frames, 1, 0.f);

	//! Mask: true for pixels that still need to be processed
	Video<char> mask(sz.width, sz.height, sz.frames, 1, false);

	//! There's a border added only if the crop
	//  doesn't touch the source image boundary
	bool border_x0 = p_crop.origin_x > 0;
	bool border_y0 = p_crop.origin_y > 0;
	bool border_t0 = p_crop.origin_t > 0;
	bool border_x1 = p_crop.ending_x < p_crop.source_sz.width;
	bool border_y1 = p_crop.ending_y < p_crop.source_sz.height;
	bool border_t1 = p_crop.ending_t < p_crop.source_sz.frames;

	//! Only pixels of the center of the image must be processed (not the boundaries)
	int n_groups = 0;
	unsigned stepx = p_params.offSet;
	unsigned stepy = p_params.offSet;
	unsigned stepf = p_params.offSetTime;
	int ori_x =                        border_x0 ? sPx-1 + sWx/2 : 0 ;
	int ori_y =                        border_y0 ? sPx-1 + sWx/2 : 0 ;
	int ori_f =                        border_t0 ? sPt-1 + sWt/2 : 0 ;
	int end_x = (int)sz.width  - (int)(border_x1 ? sPx-1 + sWx/2 : sPx-1);
	int end_y = (int)sz.height - (int)(border_y1 ? sPx-1 + sWx/2 : sPx-1);
	int end_f = (int)sz.frames - (int)(border_t1 ? sPt-1 + sWt/2 : sPt-1);
#ifdef SPTWO
	ori_f =                        border_t0 ? sPt/2-1 + sWt/2 : sPt/2 ;
	end_f = (int)sz.frames - (int)(border_t1 ? sPt/2-1 + sWt/2 : sPt/2);
#endif

	if (p_params.onlyFrame >=0)
	{
		ori_f = p_params.onlyFrame;
		end_f = p_params.onlyFrame + 1;
	}

//	ori_f = std::max((int)sz.frames / 2 - (int)sPt - 1, 0);
//	end_f = sz.frames / 2 + 1;

	for (int f = ori_f, df = 0; f < end_f; f++, df++)
	for (int y = ori_y, dy = 0; y < end_y; y++, dy++)
	for (int x = ori_x, dx = 0; x < end_x; x++, dx++)
	{
		if ( (df % stepf == 0) || (!border_t1 && f == end_f - 1))
		{
			int phasey = (!border_t1 && f == end_f - 1) ? 0 : f/stepf;

			if ( (dy % stepy == phasey % stepy) ||
			     (!border_y1 && y == end_y - 1) ||
				  (!border_y0 && y == ori_y    ) )
			{
				int phasex = (!border_y1 && y == end_y - 1) ? 0 : (phasey + y/stepy);

				if ( (dx % stepx == phasex % stepx) ||
				     (!border_x1 && x == end_x - 1) ||
				     (!border_x0 && x == ori_x    ) )
				{
					mask(x,y,f) = true;
					n_groups++;
				}
			}
		}
	}

//	printf("Processing at most %d groups of similar patches\n", n_groups);


#ifdef DEBUG_SHOW_WEIGHT
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		Video<float> mask_f(mask.sz);
		for (int i = 0; i < mask.sz.whcf; ++i) mask_f(i) = 255*(float)mask(i);
		sprintf(name, "dump/msk_step%d_%d.%d.%d_%%03d.png", step1 ? 1 : 2, part_x, part_y, part_t);
		mask_f.saveVideo(name, 1, 1);
	}
#endif

	//! Used matrices during Bayes' estimate
	const unsigned patch_dim = step1 ? sPx * sPx * sPt : sPx * sPx * sPt * sz.channels;
#ifndef SPTWO
	const unsigned patch_num = sWx * sWx * sWt;
#else
	const unsigned patch_num = sWx * sWx * sWt * sPt;
#endif

	//! Matrices used for Bayes' estimate
	vector<unsigned> index(patch_num);
	matWorkspace mat;
	mat.groupTranspose.resize(patch_num * patch_dim);
	mat.tmpMat          .resize(patch_dim * patch_dim);
	mat.covMat          .resize(patch_dim * patch_dim);
	mat.covMatTmp       .resize(patch_dim * patch_dim);
	mat.baricenter      .resize(patch_dim);

	//! Variance captured by the principal components
	Video<float> variance(mask.sz);

	//! Total number of groups of similar patches processed
	unsigned group_counter = 0;

#ifdef DEBUG_SHOW_GROUPS
	//! open group files for writing
	char n[1024];
	sprintf(n, "groups-nisy%d.dat", step1?1:2); FILE *file_groups_nisy = fopen(n,"w");
	sprintf(n, "groups-deno%d.dat", step1?1:2); FILE *file_groups_deno = fopen(n,"w");
#endif

	if (step1)
	{
		//! Allocate Sizes
		io_imBasic.resize(sz);

		//! Matrices used for Bayes' estimate
		vector<vector<float> > group(sz.channels, vector<float>(patch_num * patch_dim));

		int remaining_groups = n_groups;
		for (unsigned pt = 0; pt < sz.frames; pt++)
		for (unsigned py = 0; py < sz.height; py++)
		for (unsigned px = 0; px < sz.width ; px++)
			if (mask(px,py,pt)) //< Only non-seen patches are processed
			{
				group_counter++;

				const unsigned ij  = sz.index(px,py,pt);
				const unsigned ij3 = sz.index(px,py,pt, 0);
				//const unsigned ij3 = (ij / sz.wh) * sz.whc + ij % sz.wh;

				if (p_params.verbose && (group_counter % 100 == 0))
				{
					int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
					int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
					               p_crop.tile_y * p_crop.ntiles_x + 
					               p_crop.tile_x;

					printf("\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE", ntiles - part_idx,
							p_crop.tile_x, p_crop.tile_y, p_crop.tile_t,
							100.f - (float)remaining_groups/(float)(n_groups)*100.f,
							ntiles - part_idx);

					std::cout << std::flush;
				}

				//! Search for similar patches around the reference one
				unsigned nSimP = estimateSimilarPatchesStep1(i_imNoisy, i_fflow, i_bflow,
						group, index, ij3, p_params, i_imClean);

				//! If we use the homogeneous area trick
				bool flatPatch = false;
				if (p_params.useHomogeneousArea)
					flatPatch = computeHomogeneousAreaStep1(group, sPx * sPx * sPt,
							patch_num, threshold, sz);

#ifdef DEBUG_SHOW_GROUPS
#ifdef FAT_ORIGINAL
				if (flatPatch == false)
#endif
				{
					//! print number of patches and patch dimension
					fprintf(file_groups_nisy, "%d %% n\n", nSimP);
					fprintf(file_groups_nisy, "%d %% d\n", patch_dim*sz.channels);
					//! print groups
					for(int k = 0; k < nSimP; k++)
					{
						//! patch coordinates
						const unsigned ind = index[k];
						unsigned px0, py0, pt0, pc0;
						sz.coords(ind, px0, py0, pt0, pc0);

						fprintf(file_groups_nisy, "%d %% px\n", px0);
						fprintf(file_groups_nisy, "%d %% py\n", py0);
						fprintf(file_groups_nisy, "%d %% pt\n", pt0);
						for(int c = 0; c < sz.channels; ++c)
						for(int i = 0; i < sPx*sPx*sPt; ++i)
							fprintf(file_groups_nisy, "%.16g \n", group[c][i*nSimP + k]);
					}
				}
#endif

				//! Bayesian estimate
#ifdef FAT_ORIGINAL
				// The Bayesian estimation is skipped with the original Flat Area
				// trick, since the denoising has been done already in the 
				// computeHomogeneousArea function.
				if (flatPatch == false)
#endif
				float variance = 200.f * (p_params.rank <= patch_dim)
					? computeBayesEstimateStep1_LR(group, mat, nInverseFailed,
							p_params, nSimP)
					: computeBayesEstimateStep1_FR(group, mat, nInverseFailed,
							p_params, nSimP);

#ifdef DEBUG_SHOW_GROUPS
#ifdef FAT_ORIGINAL
				if (flatPatch == false)
#endif
				{
					//! print number of patches and patch dimension
					fprintf(file_groups_deno, "%d %% n\n", nSimP);
					fprintf(file_groups_deno, "%d %% d\n", patch_dim*sz.channels);
					//! print groups
					for(int k = 0; k < nSimP; k++)
					{
						//! patch coordinates
						const unsigned ind = index[k];
						unsigned px0, py0, pt0, pc0;
						sz.coords(ind, px0, py0, pt0, pc0);

						fprintf(file_groups_deno, "%d %% px\n", px0);
						fprintf(file_groups_deno, "%d %% py\n", py0);
						fprintf(file_groups_deno, "%d %% pt\n", pt0);
						for(int c = 0; c < sz.channels; ++c)
						for(int i = 0; i < sPx*sPx*sPt; ++i)
							fprintf(file_groups_deno, "%.16g \n", group[c][i*nSimP + k]);
					}
				}
#endif

				//! Aggregation
				remaining_groups -=
					computeAggregationStep1(io_imBasic, weight, mask, group, index,
						p_params, nSimP);
			}

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, io_imBasic, weight);

		if (p_params.verbose)
		{
			int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
			int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
			               p_crop.tile_y * p_crop.ntiles_x +
			               p_crop.tile_x;

			printf(ANSI_GRN "\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE" ANSI_RST,
					ntiles - part_idx, p_crop.tile_x, p_crop.tile_y,
					p_crop.tile_t, 100.f, ntiles - part_idx);

			std::cout << std::flush;

			//printf("Processed %d groups\n", group_counter);
		}
	}
	else
	{
		//! Allocate Sizes
		o_imFinal.resize(sz);

		//! Matrices used for Bayes' estimate
		vector<float> groupNoisy(patch_num * patch_dim);
		vector<float> groupBasic(patch_num * patch_dim);

		int remaining_groups = n_groups;
		for (unsigned pt = 0; pt < sz.frames; pt++)
		for (unsigned py = 0; py < sz.height; py++)
		for (unsigned px = 0; px < sz.width ; px++)
			if (mask(px,py,pt)) //< Only non-seen patches are processed
			{
				group_counter++;

				const unsigned ij  = sz.index(px,py,pt);
				const unsigned ij3 = sz.index(px,py,pt, 0);
				//const unsigned ij3 = (ij / sz.wh) * sz.whc + ij % sz.wh;

				if (p_params.verbose && (group_counter % 100 == 0))
				{
					int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
					int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
					               p_crop.tile_y * p_crop.ntiles_x +
					               p_crop.tile_x;

					printf("\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE", ntiles - part_idx,
							p_crop.tile_x, p_crop.tile_y, p_crop.tile_t,
							100.f - (float)remaining_groups/(float)(n_groups)*100.f,
							ntiles - part_idx);

					std::cout << std::flush;
				}

				//! Search for similar patches around the reference one
				unsigned nSimP = estimateSimilarPatchesStep2(i_imNoisy, io_imBasic,
				      i_fflow, i_bflow, groupNoisy, groupBasic, index, ij3, p_params, 
				      i_imClean);

				//! If we use the homogeneous area trick
				bool flatPatch = false;
				if (p_params.useHomogeneousArea)
					flatPatch = computeHomogeneousAreaStep2(groupNoisy,
							groupBasic, sPx * sPx *sPt, nSimP, threshold, sz);

#ifdef DEBUG_SHOW_GROUPS
#ifdef FAT_ORIGINAL
				if (flatPatch == false)
#endif
				{
					//! print number of patches and patch dimension
					fprintf(file_groups_nisy, "%d %% n\n", nSimP);
					fprintf(file_groups_nisy, "%d %% d\n", patch_dim);
					//! print groups
					for(int k = 0; k < nSimP; k++)
					{
						//! patch coordinates
						const unsigned ind = index[k];
						unsigned px0, py0, pt0, pc0;
						sz.coords(ind, px0, py0, pt0, pc0);

						fprintf(file_groups_nisy, "%d %% px\n", px0);
						fprintf(file_groups_nisy, "%d %% py\n", py0);
						fprintf(file_groups_nisy, "%d %% pt\n", pt0);
						for(int i = 0; i < patch_dim; ++i)
							fprintf(file_groups_nisy, "%.16g \n", groupNoisy[i*nSimP + k]);
					}
				}
#endif

				//! Bayesian estimate
#ifdef FAT_ORIGINAL
				// The Bayesian estimation is skipped with the original Flat Area
				// trick, since the denoising has been done already in the 
				// computeHomogeneousArea function.
				if (flatPatch == false)
#endif
				variance(ij) = 200.f * (p_params.rank <= patch_dim)
					? computeBayesEstimateStep2_LR(groupNoisy, groupBasic, mat,
							nInverseFailed, sz, p_params, nSimP, flatPatch)
					: computeBayesEstimateStep2_FR(groupNoisy, groupBasic, mat,
							nInverseFailed, sz, p_params, nSimP);

#ifdef DEBUG_SHOW_GROUPS
#ifdef FAT_ORIGINAL
				if (flatPatch == false)
#endif
				{
					//! print number of patches and patch dimension
					fprintf(file_groups_deno, "%d %% n\n", nSimP);
					fprintf(file_groups_deno, "%d %% d\n", patch_dim);
					//! print groups
					for(int k = 0; k < nSimP; k++)
					{
						//! patch coordinates
						const unsigned ind = index[k];
						unsigned px0, py0, pt0, pc0;
						sz.coords(ind, px0, py0, pt0, pc0);

						fprintf(file_groups_deno, "%d %% px\n", px0);
						fprintf(file_groups_deno, "%d %% py\n", py0);
						fprintf(file_groups_deno, "%d %% pt\n", pt0);
						for(int i = 0; i < patch_dim; ++i)
							fprintf(file_groups_deno, "%.16g \n", groupNoisy[i*nSimP + k]);
					}
				}
#endif

				//! Aggregation
				remaining_groups -=
					computeAggregationStep2(o_imFinal, weight, mask, groupNoisy,
						variance, index, p_params, nSimP);
			}

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, o_imFinal, weight);

		if (p_params.verbose)
		{
			int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
			int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
			               p_crop.tile_y * p_crop.ntiles_x +
			               p_crop.tile_x;

			printf(ANSI_GRN "\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE" ANSI_RST,
					ntiles - part_idx, p_crop.tile_x, p_crop.tile_y,
					p_crop.tile_t, 100.f, ntiles - part_idx);

			std::cout << std::flush;
		}
	}


	if (nInverseFailed > 0 && p_params.verbose)
		std::cout << "nInverseFailed = " << nInverseFailed << std::endl;

#ifdef DEBUG_SHOW_WEIGHT
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		sprintf(name, "dump/weight_step%d_%d.%d.%d_%%03d.png",
				step1 ? 1 : 2, part_x, part_y, part_t);
		weight.saveVideo(name, 1, 1);
	}
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		sprintf(name, "dump/var_step%d_%d.%d.%d_%%03d.png",
				step1 ? 1 : 2, part_x, part_y, part_t);
		variance.saveVideo(name, 1, 1);
	}
#endif

#ifdef DEBUG_SHOW_GROUPS
	//! Close group files
	fclose(file_groups_nisy);
	fclose(file_groups_deno);
#endif

	return group_counter;
}

#ifndef SPTWO
/**
 * @brief Keep from all near patches the similar ones to the reference patch
 * for the first step. The input optical flow is used to center the search
 *
 * windows along the trajectory of the reference patch.
 *
 * @param i_im: contains the noisy image on which distances are processed;
 * @param i_fflow: forward optical flow (optional)
 * @param i_bflow: backward optical flow (optional)
 * @param o_group: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
unsigned estimateSimilarPatchesStep1(
	Video<float> const& im
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<std::vector<float> > &group
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
,  Video<float> const &imClean
){
	//! Initialization
	int sWx   = params.sizeSearchWindow;
	int sWy   = params.sizeSearchWindow;
	const int sWt_f = params.sizeSearchTimeRangeFwd;
	const int sWt_b = params.sizeSearchTimeRangeBwd;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	const VideoSize sz = im.sz;
	const bool use_flow = (fflow.sz.width > 0);

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	sz.coords(pidx, px, py, pt, pc);

	//! Temporal search range
	int ranget[2];
#ifdef CENTRED_SEARCH
	ranget[0] = std::max(0, (int)pt -  sWt_b);
	ranget[1] = std::min((int)sz.frames - sPt, (int)pt + sWt_f);
#else
	int shift_t = std::min(0, (int)pt -  sWt_b)  
	            + std::max(0, (int)pt +  sWt_f - (int)sz.frames + sPt); 

	ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
	ranget[1] = std::min((int)sz.frames - sPt, (int)pt +  sWt_f - shift_t);
#endif

	//! Redefine size of temporal search range
	int sWt = ranget[1] - ranget[0] + 1;

	//! Allocate vector of patch distances
	std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

	//! Number of patches in search region
	int nsrch = 0;

	//! Schedule frames in temporal range
	std::vector<int> srch_ranget;
	srch_ranget.push_back(pt);
	for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
	for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

	//! Trajectory of search center
	std::vector<int> cx(sWt,0), cy(sWt,0), ct(sWt,0);

	//! Search
	for (int ii = 0; ii < sWt; ++ii)
	{
		int qt = srch_ranget[ii]; // video frame number
		int dt = qt - ranget[0]; // search region frame number
		int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

		//! Integrate optical flow to new center
		if (dir != 0)
		{
			int cx0 = cx[dt - dir];
			int cy0 = cy[dt - dir];
			int ct0 = ct[dt - dir];
			float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
			float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
			cx[dt] = std::max(0., std::min((double)sz.width  - 1, round(cx_f)));
			cy[dt] = std::max(0., std::min((double)sz.height - 1, round(cy_f)));
			ct[dt] = qt;
		}
		else
		{
			cx[dt] = px;
			cy[dt] = py;
			ct[dt] = pt;
		}

		//! Spatial search range
		int rangex[2];
		int rangey[2];
#ifdef CENTRED_SEARCH
		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2);
#else
		int shift_x = std::min(0, cx[dt] - (sWx-1)/2); 
		int shift_y = std::min(0, cy[dt] - (sWy-1)/2); 

		shift_x += std::max(0, cx[dt] + (sWx-1)/2 - (int)sz.width  + sPx); 
		shift_y += std::max(0, cy[dt] + (sWy-1)/2 - (int)sz.height + sPx); 

		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2 - shift_x);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2 - shift_y);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2 - shift_x);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2 - shift_y);
#endif

		//! Check if oracle has been provided
		const Video<float> *p_im = (imClean.sz.whcf) ? &imClean : &im;

		//! Compute distance between patches in search range
		for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance
			float dist = 0.f, dif;
			for (int ht = 0; ht < sPt; ht++)
			for (int hy = 0; hy < sPx; hy++)
			for (int hx = 0; hx < sPx; hx++)
				dist += (dif = (*p_im)(px + hx, py + hy, pt + ht)
				             - (*p_im)(qx + hx, qy + hy, qt + ht)) * dif;
	
			//! Save distance and corresponding patch index
			distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
		}
	}

	distance.resize(nsrch);

	//! Keep only the N2 best similar patches
	unsigned nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nSimP,
	                  distance.end(), comparaisonFirst);

	if (nSimP <  params.nSimilarPatches)
		printf("SR1 [%d,%d,%d] ~ nsim = %d\n", px,py,pt,nSimP);

	//! Add more patches if their distance is 0
	// FIXME: use a distance threshold as in the 2nd step
	const float threshold = distance[nSimP - 1].first;
	nSimP = 0;
	for (unsigned n = 0; n < distance.size(); n++)
		if (distance[n].first <= threshold)
			index[nSimP++] = distance[n].second;
	for (int n = 0; n < nSimP; n++) index[n] = distance[n].second;

	//! Stack selected patches into the 3D group
	const int w   = sz.width;
	const int wh  = sz.wh;
	const int whc = sz.whc;
	for (int c  = 0; c < sz.channels; c++)
	for (int ht = 0, k = 0; ht < sPt; ht++)
	for (int hy = 0;        hy < sPx; hy++)
	for (int hx = 0;        hx < sPx; hx++)
	for (int n  = 0; n < nSimP; n++, k++)
		group[c][k] = im(c * wh + index[n] + ht * whc + hy * w + hx);

	/* 000  pixels from all patches
	 * 001  pixels from all patches
	 * ...
	 * spt,spx,spx pixels from all patches
	 */

	return nSimP;
}
#else // SPTWO
/**
 * @brief Keep from all near patches the similar ones to the reference patch
 * for the first step. The input optical flow is used to center the search
 * windows along the trajectory of the reference patch.
 *
 * @param i_im: contains the noisy image on which distances are processed;
 * @param i_fflow: forward optical flow (optional)
 * @param i_bflow: backward optical flow (optional)
 * @param o_group: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
unsigned estimateSimilarPatchesStep1(
	Video<float> const& im
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<std::vector<float> > &group
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
,  Video<float> const &imClean
){
	//! Initialization
	int sWx   = params.sizeSearchWindow;
	int sWy   = params.sizeSearchWindow;
	const int sWt_f = params.sizeSearchTimeRangeFwd;
	const int sWt_b = params.sizeSearchTimeRangeBwd;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime/2;
	const VideoSize sz = im.sz;
	const bool use_flow = (fflow.sz.width > 0);

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	sz.coords(pidx, px, py, pt, pc);

	//! Temporal search range
	int ranget[2];
#ifdef CENTRED_SEARCH
	ranget[0] = std::max(sPt, (int)pt - sWt_b);
	ranget[1] = std::min((int)sz.frames - sPt - 1, (int)pt + sWt_f);
#else
	int shift_t = std::min(0, (int)pt - sWt_b - sPt)  
	            + std::max(0, (int)pt + sWt_f - (int)sz.frames + sPt + 1); 

	ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
	ranget[1] = std::min((int)sz.frames - sPt - 1, (int)pt +  sWt_f - shift_t);
#endif

	//! Redefine size of temporal search range
	int sWt = ranget[1] - ranget[0] + 1;

	//! Allocate vector of patch distances
	std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

	//! Number of patches in search region
	int nsrch = 0;

	//! Schedule frames in temporal range
	std::vector<int> srch_ranget;
	srch_ranget.push_back(pt);
	for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
	for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

	//! Trajectory of search center
	std::vector<int> cx(sWt,0), cy(sWt,0), ct(sWt,0);

	//! Search
	for (int ii = 0; ii < sWt; ++ii)
	{
		int qt = srch_ranget[ii]; // video frame number
		int dt = qt - ranget[0]; // search region frame number
		int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

		//! Integrate optical flow to new center
		if (dir != 0)
		{
			int cx0 = cx[dt - dir];
			int cy0 = cy[dt - dir];
			int ct0 = ct[dt - dir];
			float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
			float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
			cx[dt] = std::max(0., std::min((double)sz.width  - 1, round(cx_f)));
			cy[dt] = std::max(0., std::min((double)sz.height - 1, round(cy_f)));
			ct[dt] = qt;
		}
		else
		{
			cx[dt] = px;
			cy[dt] = py;
			ct[dt] = pt;
		}

		//! Spatial search range
		int rangex[2];
		int rangey[2];
#ifdef CENTRED_SEARCH
		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2);
#else
		int shift_x = std::min(0, cx[dt] - (sWx-1)/2); 
		int shift_y = std::min(0, cy[dt] - (sWy-1)/2); 

		shift_x += std::max(0, cx[dt] + (sWx-1)/2 - (int)sz.width  + sPx); 
		shift_y += std::max(0, cy[dt] + (sWy-1)/2 - (int)sz.height + sPx); 

		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2 - shift_x);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2 - shift_y);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2 - shift_x);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2 - shift_y);
#endif

		//! Check if oracle has been provided
		const Video<float> *p_im = (imClean.sz.whcf) ? &imClean : &im;

		//! Compute distance between patches in search range
		for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance
			float dist = 0.f, dif;

			//! Forward patch trajectory
			int t_px = px, t_qx = qx;
			int t_py = py, t_qy = qy;
			int t_pt = pt, t_qt = qt;
			for (int ht = 0; ht <= sPt; ht++)
			{
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = (*p_im)(BC_SYMM, t_px + hx, t_py + hy, t_pt)
					             - (*p_im)(BC_SYMM, t_qx + hx, t_qy + hy, t_qt)) * dif;

				t_px = t_px + (use_flow ? round(fflow(t_px,t_py,t_pt,0)) : 0);
				t_qx = t_qx + (use_flow ? round(fflow(t_qx,t_qy,t_qt,0)) : 0);
				t_py = t_py + (use_flow ? round(fflow(t_px,t_py,t_pt,1)) : 0);
				t_qy = t_qy + (use_flow ? round(fflow(t_qx,t_qy,t_qt,1)) : 0);

				t_px = std::max(0, std::min(t_px, (int)sz.width  - 1));
				t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
				t_py = std::max(0, std::min(t_py, (int)sz.height - 1));
				t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));

				t_pt++;
				t_qt++;
			}

			//! Backward patch trajectory
			t_px = px, t_qx = qx;
			t_py = py, t_qy = qy;
			t_pt = pt, t_qt = qt;
			for (int ht = -1; ht >= -sPt; ht--)
			{
				t_px = t_px + (use_flow ? round(bflow(t_px,t_py,t_pt,0)) : 0);
				t_qx = t_qx + (use_flow ? round(bflow(t_qx,t_qy,t_qt,0)) : 0);
				t_py = t_py + (use_flow ? round(bflow(t_px,t_py,t_pt,1)) : 0);
				t_qy = t_qy + (use_flow ? round(bflow(t_qx,t_qy,t_qt,1)) : 0);

				t_px = std::max(0, std::min(t_px, (int)sz.width  - 1));
				t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
				t_py = std::max(0, std::min(t_py, (int)sz.height - 1));
				t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));

				t_pt--;
				t_qt--;

				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = im(BC_SYMM, t_px + hx, t_py + hy, t_pt)
									 - im(BC_SYMM, t_qx + hx, t_qy + hy, t_qt)) * dif;
			}
	
			//! Save distance and corresponding patch index
			distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
		}
	}

	distance.resize(nsrch);

	//! Keep only the n2/(2*sPt + 1) most similar extended patches
	unsigned nSimP = std::min(params.nSimilarPatches/(2*sPt + 1), (unsigned)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nSimP,
	                  distance.end(), comparaisonFirst);

	if (nSimP < params.nSimilarPatches/(2*sPt + 1))
		printf("SR1 [%d,%d,%d] ~ nsim = %d\n", px,py,pt,nSimP);

	//! Register position of 2D patches
	index.resize(nSimP * (2*sPt + 1));
	for (int n = 0, nn = 0; n < nSimP; n++)
	{
		unsigned t_qx0, t_qy0, t_qt0, t_qc0;
		sz.coords(distance[n].second, t_qx0, t_qy0, t_qt0, t_qc0);

		//! Forward patch trajectory
		int t_qx = t_qx0;
		int t_qy = t_qy0;
		int t_qt = t_qt0;
		for (int ht = 0; ht <= sPt; ht++)
		{
/*			if(t_qt0 == 20)
			{
				printf("\tt_q = (%d,%d,%d)   ", t_qx, t_qy, t_qt);
				printf("\tf(t_q) = (%f,%f) \n", fflow(t_qx,t_qy,t_qt,0), fflow(t_qx,t_qy,t_qt,1));
			}*/

			index[nn++] = sz.index(t_qx, t_qy, t_qt, 0);
			t_qx = t_qx + (use_flow ? round(fflow(t_qx,t_qy,t_qt,0)) : 0);
			t_qy = t_qy + (use_flow ? round(fflow(t_qx,t_qy,t_qt,1)) : 0);
			t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
			t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));
			t_qt++;
		}

		//! Backward patch trajectory
		t_qx = t_qx0;
		t_qy = t_qy0;
		t_qt = t_qt0;
		for (int ht = -1; ht >= -sPt; ht--)
		{
			t_qx = t_qx + (use_flow ? round(bflow(t_qx,t_qy,t_qt,0)) : 0);
			t_qy = t_qy + (use_flow ? round(bflow(t_qx,t_qy,t_qt,1)) : 0);
			t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
			t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));
			t_qt--;
			index[nn++] = sz.index(t_qx, t_qy, t_qt);

/*			if(t_qt0 == 70)
			{
				printf("\tt_q = (%d,%d,%d)   ", t_qx, t_qy, t_qt);
				printf("\tf(t_q) = (%f,%f) \n", bflow(t_qx,t_qy,t_qt,0), bflow(t_qx,t_qy,t_qt,1));
			}*/

		}
	}

	//! Update number of similar patches as number of 2D slices
	nSimP *= 2*sPt + 1; 

	//! Stack selected 2D patches into the 3D group
	const int w   = sz.width;
	const int wh  = sz.wh;
	const int whc = sz.whc;
	unsigned qx, qy, qt, qc;
	for (int c  = 0; c < sz.channels; c++)
	for (int hy = 0, k = 0; hy < sPx; hy++)
	for (int hx = 0;        hx < sPx; hx++)
//	for (int ht = 0;        ht < sPt; ht++)
	for (int n  = 0; n < nSimP; n++, k++)
	{
		sz.coords(index[n], qx, qy, qt, qc);
		group[c][k] = im(BC_SYMM, qx + hx, qy + hy, qt, qc + c);
	}

	/* 000  pixels from all patches
	 * 001  pixels from all patches
	 * ...
	 * spt,spx,spx pixels from all patches
	 */

	return nSimP;
}
#endif

#ifndef SPTWO
/**
 * @brief Keep from all near patches the similar ones to the reference patch
 * for the second step. The input optical flow is used to center the search
 * windows along the trajectory of the reference patch.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param i_fflow: forward optical flow (optional)
 * @param i_bflow: backward optical flow (optional)
 * @param o_groupNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_groupBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of images;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
,  Video<float> const &imClean
){
	//! Initialization
	int sWx   = params.sizeSearchWindow;
	int sWy   = params.sizeSearchWindow;
	const int sWt_f = params.sizeSearchTimeRangeFwd;
	const int sWt_b = params.sizeSearchTimeRangeBwd;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	const VideoSize sz = imNoisy.sz;
	const int chnls = sz.channels;
	bool use_flow = (fflow.sz.width > 0);

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	sz.coords(pidx, px, py, pt, pc);

	//! Temporal search range
	int ranget[2];
#ifdef CENTRED_SEARCH
	ranget[0] = std::max(0, (int)pt -  sWt_b);
	ranget[1] = std::min((int)sz.frames - sPt, (int)pt + sWt_f);
#else
	int shift_t = std::min(0, (int)pt -  sWt_b)  
	            + std::max(0, (int)pt +  sWt_f - (int)sz.frames + sPt); 

	ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
	ranget[1] = std::min((int)sz.frames - sPt, (int)pt +  sWt_f - shift_t);
#endif

	//! Redefine size of temporal search range
	int sWt = ranget[1] - ranget[0] + 1;

	//! Allocate vector of patch distances
	std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

	//! Number of patches in search region
	int nsrch = 0;

	//! Schedule frames in temporal range
	std::vector<int> srch_ranget;
	srch_ranget.push_back(pt);
	for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
	for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

	//! Trajectory of search center
	std::vector<int> cx(sWt,0), cy(sWt,0), ct(sWt,0);

	//! Search
	for (int ii = 0; ii < sWt; ++ii)
	{
		int qt = srch_ranget[ii]; // video frame number
		int dt = qt - ranget[0]; // search region frame number
		int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

		//! Integrate optical flow to new center
		if (dir != 0)
		{
			int cx0 = cx[dt - dir];
			int cy0 = cy[dt - dir];
			int ct0 = ct[dt - dir];
//			if (cx0 > (int)bflow.sz.width - 1 || cx0 < 0 ||
//			    cy0 > (int)bflow.sz.height- 1 || cy0 < 0 ||
//			    ct0 > (int)bflow.sz.frames- 1 || ct0 < 0 )
//				printf("sz = %d,%d,%d ~ c = %d,%d,%d\n",bflow.sz.width, bflow.sz.height, bflow.sz.frames, cx0, cy0, ct0);
			float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
			float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
			cx[dt] = std::max(0., std::min((double)sz.width  - 1, round(cx_f)));
			cy[dt] = std::max(0., std::min((double)sz.height - 1, round(cy_f)));
			ct[dt] = qt;
		}
		else
		{
			cx[dt] = px;
			cy[dt] = py;
			ct[dt] = pt;
		}

		//! Spatial search range
		int rangex[2];
		int rangey[2];
#ifdef CENTRED_SEARCH
		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2);
#else
		int shift_x = std::min(0, cx[dt] - (sWx-1)/2); 
		int shift_y = std::min(0, cy[dt] - (sWy-1)/2); 

		shift_x += std::max(0, cx[dt] + (sWx-1)/2 - (int)sz.width  + sPx); 
		shift_y += std::max(0, cy[dt] + (sWy-1)/2 - (int)sz.height + sPx); 

		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2 - shift_x);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2 - shift_y);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2 - shift_x);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2 - shift_y);
#endif

		//! Check if oracle has been provided
		const Video<float> *p_im = (imClean.sz.whcf) ? &imClean : &imBasic;

		//! Compute distance between patches in search range
		for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance
			float dist = 0.f, dif;
			for (int c = 0; c < chnls; c++)
			for (int ht = 0; ht < sPt; ht++)
			for (int hy = 0; hy < sPx; hy++)
			for (int hx = 0; hx < sPx; hx++)
				dist += (dif = (*p_im)(px + hx, py + hy, pt + ht, c)
				             - (*p_im)(qx + hx, qy + hy, qt + ht, c) ) * dif;
	
			//! Save distance and corresponding patch index
			distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
		}
	}

	distance.resize(nsrch);

	//! Keep only the nSimilarPatches best similar patches
	unsigned nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nSimP,
	                  distance.end(), comparaisonFirst);

	if (nSimP <  params.nSimilarPatches)
		printf("SR2 [%d,%d,%d] ~ nsim = %d\n", px,py,pt,nSimP);

	//! Add more patches if their distance is bellow the threshold
	const float threshold = (params.tau > distance[nSimP - 1].first ?
	                         params.tau : distance[nSimP - 1].first);
	nSimP = 0;
	for (unsigned n = 0; n < distance.size(); n++)
		if (distance[n].first <= threshold)
			index[nSimP++] = distance[n].second;

	//! Save similar patches into 3D groups
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;
	for (unsigned c = 0, k = 0; c < chnls; c++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++)
	for (unsigned n = 0; n < nSimP; n++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
	}

	return nSimP;
}
#else // SPTWO
/**
 * @brief Keep from all near patches the similar ones to the reference patch
 * for the second step. The input optical flow is used to center the search
 * windows along the trajectory of the reference patch.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param i_fflow: forward optical flow (optional)
 * @param i_bflow: backward optical flow (optional)
 * @param o_groupNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_groupBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of images;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
,  Video<float> const &imClean
){
	//! Initialization
	int sWx   = params.sizeSearchWindow;
	int sWy   = params.sizeSearchWindow;
	const int sWt_f = params.sizeSearchTimeRangeFwd;
	const int sWt_b = params.sizeSearchTimeRangeBwd;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime/2;
	const VideoSize sz = imNoisy.sz;
	const int chnls = sz.channels;
	bool use_flow = (fflow.sz.width > 0);

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	sz.coords(pidx, px, py, pt, pc);

	//! Temporal search range
	int ranget[2];
#ifdef CENTRED_SEARCH
	ranget[0] = std::max(sPt, (int)pt - sWt_b);
	ranget[1] = std::min((int)sz.frames - sPt - 1, (int)pt + sWt_f);
#else
	int shift_t = std::min(0, (int)pt - sWt_b - sPt)  
	            + std::max(0, (int)pt + sWt_f - (int)sz.frames + sPt + 1); 

	ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
	ranget[1] = std::min((int)sz.frames - sPt - 1, (int)pt +  sWt_f - shift_t);
#endif

	//! Redefine size of temporal search range
	int sWt = ranget[1] - ranget[0] + 1;

	//! Allocate vector of patch distances
	std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

	//! Number of patches in search region
	int nsrch = 0;

	//! Schedule frames in temporal range
	std::vector<int> srch_ranget;
	srch_ranget.push_back(pt);
	for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
	for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

	//! Trajectory of search center
	std::vector<int> cx(sWt,0), cy(sWt,0), ct(sWt,0);

	//! Search
	for (int ii = 0; ii < sWt; ++ii)
	{
		int qt = srch_ranget[ii]; // video frame number
		int dt = qt - ranget[0]; // search region frame number
		int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

		//! Integrate optical flow to new center
		if (dir != 0)
		{
			int cx0 = cx[dt - dir];
			int cy0 = cy[dt - dir];
			int ct0 = ct[dt - dir];
//			if (cx0 > (int)bflow.sz.width - 1 || cx0 < 0 ||
//			    cy0 > (int)bflow.sz.height- 1 || cy0 < 0 ||
//			    ct0 > (int)bflow.sz.frames- 1 || ct0 < 0 )
//				printf("sz = %d,%d,%d ~ c = %d,%d,%d\n",bflow.sz.width, bflow.sz.height, bflow.sz.frames, cx0, cy0, ct0);
			float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
			float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
			cx[dt] = std::max(0., std::min((double)sz.width  - 1, round(cx_f)));
			cy[dt] = std::max(0., std::min((double)sz.height - 1, round(cy_f)));
			ct[dt] = qt;
		}
		else
		{
			cx[dt] = px;
			cy[dt] = py;
			ct[dt] = pt;
		}

		//! Spatial search range
		int rangex[2];
		int rangey[2];
#ifdef CENTRED_SEARCH
		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2);
#else
		int shift_x = std::min(0, cx[dt] - (sWx-1)/2); 
		int shift_y = std::min(0, cy[dt] - (sWy-1)/2); 

		shift_x += std::max(0, cx[dt] + (sWx-1)/2 - (int)sz.width  + sPx); 
		shift_y += std::max(0, cy[dt] + (sWy-1)/2 - (int)sz.height + sPx); 

		rangex[0] = std::max(0, cx[dt] - (sWx-1)/2 - shift_x);
		rangey[0] = std::max(0, cy[dt] - (sWy-1)/2 - shift_y);

		rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2 - shift_x);
		rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2 - shift_y);
#endif

		//! Check if oracle has been provided
		const Video<float> *p_im = (imClean.sz.whcf) ? &imClean : &imBasic;

		//! Compute distance between patches in search range
		for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance
			float dist = 0.f, dif;

			//! Forward patch trajectory
			int t_px = px, t_qx = qx;
			int t_py = py, t_qy = qy;
			int t_pt = pt, t_qt = qt;
			for (int ht = 0; ht <= sPt; ht++)
			{
				for (int c  = 0; c < chnls; c++)
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = (*p_im)(BC_SYMM, t_px + hx, t_py + hy, t_pt, c)
					             - (*p_im)(BC_SYMM, t_qx + hx, t_qy + hy, t_qt, c)) * dif;

				t_px = t_px + (use_flow ? round(fflow(t_px,t_py,t_pt,0)) : 0);
				t_qx = t_qx + (use_flow ? round(fflow(t_qx,t_qy,t_qt,0)) : 0);
				t_py = t_py + (use_flow ? round(fflow(t_px,t_py,t_pt,1)) : 0);
				t_qy = t_qy + (use_flow ? round(fflow(t_qx,t_qy,t_qt,1)) : 0);

				t_px = std::max(0, std::min(t_px, (int)sz.width  - 1));
				t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
				t_py = std::max(0, std::min(t_py, (int)sz.height - 1));
				t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));

				t_pt++;
				t_qt++;
			}

			//! Backward patch trajectory
			t_px = px, t_qx = qx;
			t_py = py, t_qy = qy;
			t_pt = pt, t_qt = qt;
			for (int ht = -1; ht >= -sPt; ht--)
			{
				t_px = t_px + (use_flow ? round(bflow(t_px,t_py,t_pt,0)) : 0);
				t_qx = t_qx + (use_flow ? round(bflow(t_qx,t_qy,t_qt,0)) : 0);
				t_py = t_py + (use_flow ? round(bflow(t_px,t_py,t_pt,1)) : 0);
				t_qy = t_qy + (use_flow ? round(bflow(t_qx,t_qy,t_qt,1)) : 0);

				t_px = std::max(0, std::min(t_px, (int)sz.width  - 1));
				t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
				t_py = std::max(0, std::min(t_py, (int)sz.height - 1));
				t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));

				t_pt--;
				t_qt--;

				for (int c  = 0; c < chnls; c++)
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = imBasic(BC_SYMM, t_px + hx, t_py + hy, t_pt, c)
									 - imBasic(BC_SYMM, t_qx + hx, t_qy + hy, t_qt, c)) * dif;
			}

			//! Save distance and corresponding patch index
			distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
		}
	}

	distance.resize(nsrch);

	//! Keep only the nSimilarPatches best similar patches
	unsigned nSimP = std::min(params.nSimilarPatches/(2*sPt + 1), (unsigned)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nSimP,
	                  distance.end(), comparaisonFirst);

	if (nSimP <  params.nSimilarPatches/(2*sPt + 1))
		printf("SR2 [%d,%d,%d] ~ nsim = %d\n", px,py,pt,nSimP);

	//! Save index of similar patches
	const float threshold = (params.tau > distance[nSimP - 1].first ?
	                         params.tau : distance[nSimP - 1].first);
	nSimP = 0;

	//! Register position of similar patches
	for (unsigned n = 0; n < distance.size(); n++)
		if (distance[n].first < threshold)
		{
			unsigned t_qx0, t_qy0, t_qt0, t_qc0;
			sz.coords(distance[n].second, t_qx0, t_qy0, t_qt0, t_qc0);

			//! Forward patch trajectory
			int t_qx = t_qx0;
			int t_qy = t_qy0;
			int t_qt = t_qt0;
			for (int ht = 0; ht <= sPt; ht++)
			{
				index[nSimP++] = sz.index(t_qx, t_qy, t_qt, 0);
				t_qx = t_qx + (use_flow ? round(fflow(t_qx,t_qy,t_qt,0)) : 0);
				t_qy = t_qy + (use_flow ? round(fflow(t_qx,t_qy,t_qt,1)) : 0);
				t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
				t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));
				t_qt++;
			}

			//! Backward patch trajectory
			t_qx = t_qx0;
			t_qy = t_qy0;
			t_qt = t_qt0;
			for (int ht = -1; ht >= -sPt; ht--)
			{
				t_qx = t_qx + (use_flow ? round(bflow(t_qx,t_qy,t_qt,0)) : 0);
				t_qy = t_qy + (use_flow ? round(bflow(t_qx,t_qy,t_qt,1)) : 0);
				t_qx = std::max(0, std::min(t_qx, (int)sz.width  - 1));
				t_qy = std::max(0, std::min(t_qy, (int)sz.height - 1));
				t_qt--;
				index[nSimP++] = sz.index(t_qx, t_qy, t_qt);
			}
		}

	//! Save similar patches into 3D groups
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;
	unsigned qx, qy, qt, qc;
	for (unsigned c = 0, k = 0; c < chnls; c++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++)
//	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned n = 0; n < nSimP; n++, k++)
	{
		sz.coords(index[n], qx, qy, qt, qc);
		groupNoisy[k] = imNoisy(BC_SYMM, qx + hx, qy + hy, qt, qc + c);
		groupBasic[k] = imBasic(BC_SYMM, qx + hx, qy + hy, qt, qc + c);
	}

	return nSimP;
}
#endif

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group: contains for each channels values of similar patches. If an homogeneous area
 *			is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_doLinearRegression: if true, apply a linear regression to average value of pixels;
 * @param p_imSize: size of the image.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep1(
	std::vector<std::vector<float> > &io_group
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Initialization
	const unsigned N = p_sP * p_nSimP;

	//! Compute the standard deviation of the set of patches
	float variance = 0.f;
	for (unsigned c = 0; c < p_imSize.channels; c++)
		variance += computeVariance(io_group[c], p_sP, p_nSimP, 1);

	if (variance < p_threshold)
	{
#ifdef FAT_ORIGINAL
		//! If we are in an homogeneous area, estimate denoised pixels as mean
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
			double mean = 0.f;
			for (unsigned k = 0; k < N; k++) mean += io_group[c][k];
			mean /= (double) N;

			for (unsigned k = 0; k < N; k++) io_group[c][k] = mean;
		}
#endif
		return 1;
	}
	else return 0;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_groupNoisy: inputs values of similar patches for the noisy video;
 *                         if the area is classified as homogeneous, outputs the
 *                         average of all pixels in all patches.
 * @param i_groupBasic: contains values of similar patches for the basic video.
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the video.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep2(
	std::vector<float> &io_groupNoisy
,	std::vector<float> const &i_groupBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Parameters
	const unsigned sPC = p_sP * p_imSize.channels;

	//! Compute the standard deviation of the set of patches
	const float variance = computeVariance(io_groupNoisy, p_sP, p_nSimP, p_imSize.channels);

	if (variance < p_threshold)
	{
#ifdef FAT_ORIGINAL
		//! If we are in an homogeneous area, estimate denoised pixels as mean
		for (unsigned c = 0, kn = 0, kb = 0; c < p_imSize.channels; c++)
		{
			double mean = 0.f;
			for (unsigned i = 0; i < p_sP * p_nSimP; i++, kn++)
				mean += i_groupBasic[kn];
			mean /= (double)(p_sP * p_nSimP);

			for (unsigned i = 0; i < p_sP * p_nSimP; i++, kb++)
				io_groupNoisy[kb] = mean;
		}
#endif
		return 1;
	}
	else return 0;
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *    - groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
float computeBayesEstimateStep1_FR(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters
	const unsigned chnls = io_group.size();
#ifndef SPTWO
	const unsigned sP2   = p_params.sizePatch * p_params.sizePatch * p_params.sizePatchTime;
#else
	const unsigned sP2   = p_params.sizePatch * p_params.sizePatch;
#endif
	const float valDiag  = p_params.beta * p_params.sigma * p_params.sigma;

	//! Check if all patches are equal
	/* In some images might have a dark noiseless frame (for instance the image
	 * might be the result of a transformation, and a black region correspond to
	 * the pixels that where before the frame before the registration). A group
	 * of patches with 0 variance create problems afterwards. Such groups have to
	 * be avoided.
	 * */

	/* TODO
	 * This fix only addresses the problem in which all patches are equal. A more
	 * subtle situation is when some patches are equal. These patches shouldn't be
	 * denoised, since it is extremely unlekely that such an event occurs.
	 * */

	for (unsigned c = 0; c < io_group.size(); c++)
	for (unsigned j = 0; j < sP2; j++)
	{
		float v = io_group[c][j * p_nSimP];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_group[c][j * p_nSimP + i])
				goto not_equal;
	}

	//! All patches are equal ~ do nothing
	return 0.f;

not_equal:
	//! Some patches are different ~ denoise

#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Bayes estimate
	if (p_nSimP >= sP2)
	{
		for (unsigned c = 0; c < chnls; c++)
		{
			//! Center data around the baricenter
			centerData(io_group[c], i_mat.baricenter, p_nSimP, sP2);

			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(io_group[c], i_mat.covMat, p_nSimP, sP2);

			//! Bayes' Filtering
			if (inverseMatrix(i_mat.covMat, sP2) == EXIT_SUCCESS)
			{
				productMatrix(i_mat.groupTranspose, i_mat.covMat, io_group[c],
								  sP2, sP2, p_nSimP);
				for (unsigned k = 0; k < sP2 * p_nSimP; k++)
					io_group[c][k] -= valDiag * i_mat.groupTranspose[k];
			}
			else io_nInverseFailed++;

			//! Add baricenter
			for (unsigned j = 0, k = 0; j < sP2; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group[c][k] += i_mat.baricenter[j];
		}
	}
	else io_nInverseFailed++;
#else
	//! Show color-code representing patch status
	if (p_nSimP >= sP2)
	{
		float color[3];
		// green
		color[0] =  147.2243f;
		color[1] =  0;
		color[2] = -208.2066f;

		for (unsigned c = 0; c < chnls; c++)
		{
			//! Center data around the baricenter
			centerData(io_group[c], i_mat.baricenter, p_nSimP, sP2);

			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(io_group[c], i_mat.covMat, p_nSimP, sP2);

			//! Bayes' Filtering
			if (inverseMatrix(i_mat.covMat, sP2) == EXIT_FAILURE)
			{
				io_nInverseFailed++;
				// blue
				color[0] =  147.2243f;
				color[1] = -180.3122f;
				color[2] =  104.1033f;
			}
		}

		for (unsigned k = 0; k < sP2 * p_nSimP; k++)
		{
			io_group[0][k] =  color[0];
			io_group[1][k] =  color[1];
			io_group[2][k] =  color[2];
		}
	}
	else
	{
		io_nInverseFailed++;
		for (unsigned k = 0; k < sP2 * p_nSimP; k++)
		{
			// red
			io_group[0][k] = 147.2243f;
			io_group[1][k] = 180.3122f;
			io_group[2][k] = 104.1033f;
		}
	}
#endif

	return 0;
}

/**
 * @brief Implementation of computeBayesEstimateStep1_LR computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the eigenvectors/values of the data covariance matrix using LAPACK.
 *
 * See computeBayesEstimateStep1_LR for information about the arguments.
 **/
float computeBayesEstimateStep1_LR_EIG_LAPACK(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
#if defined(PAUL_VARIANCE) || defined(PAUL_SIMPLE_VARIANCE)
//	const float sigma2 = p_params.sigma * p_params.sigma;
	const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
#else
	const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
#endif
	const float thres  = p_params.variThres;
#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;
#else
	//! SPTWO consideres a Gaussian model for 2D patches
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch;
#endif
	const unsigned r    = p_params.rank;

	//! Variances
	float  rank_variance = 0.f;
	float total_variance = 0.f;

	for (unsigned c = 0; c < io_group.size(); c++)
	{
		//! Center 3D group
		centerData(io_group[c], i_mat.baricenter, p_nSimP, sPC);

		if (r > 0)
		{
			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(io_group[c], i_mat.covMat, p_nSimP, sPC);

			//! Estimate total variance
			for (int i = 0; i < sPC; ++i)
				total_variance += std::max(i_mat.covMat[i*sPC + i] - sigma2, 0.f);

			//! Compute leading eigenvectors
			int info = matrixEigs(i_mat.covMat, sPC, r, i_mat.covEigVals, i_mat.covEigVecs);

			//! Substract sigma2 and compute variance captured by the r leading eigenvectors
			int rr = 0; // final rank (rr = r, except if variThres != 0)
			for (int i = 0; i < r; ++i)
			{
#if !(defined(SPTWO) || defined(VARIANCE_THRESHOLD_1))
 #ifdef CLIPPED_VARIANCE
				/* thres is used to control a variance threshold proportional
				 * to the noise. If the eigenvalue is bellow the threshold, the 
				 * filter coefficient is set to zero. Otherwise it is set to 
				 * the clipped Wiener coefficient.
				 */
				if (thres != 0.f)
				{
					rr += (i_mat.covEigVals[i] > thres * sigma2) ? 1 : 0;
					i_mat.covEigVals[i] -= (i_mat.covEigVals[i] > thres * sigma2) ? 
					                       std::min(i_mat.covEigVals[i], sigma2) : 
					                       i_mat.covEigVals[i];
				}
				else
					i_mat.covEigVals[i] -= std::min(i_mat.covEigVals[i], sigma2);
 #elif defined(PAUL_VARIANCE)
				/* variance estimate obtained by inverting Debashis Paul asymptotic
				 * limit for the expected value of the eigenvalues.
				 */
//				float tmp, gamma = p_params.beta * (float)sPC/(float)p_nSimP;
				float tmp, gamma = (float)sPC/(float)p_nSimP;

				if (i_mat.covEigVals[i] > sigma2 * (tmp = (1 + sqrtf(gamma)))*tmp)
				{
					rr ++;
					tmp = i_mat.covEigVals[i] - sigma2 * (1 + gamma);
					i_mat.covEigVals[i] = tmp * 0.5 
						                 * (1. + sqrtf(std::max(0., 1. - 4.*gamma*sigma2*sigma2/tmp/tmp)));
				}
				else
					i_mat.covEigVals[i] = 0;

 #elif defined(PAUL_SIMPLE_VARIANCE)
				/* variance estimate obtained by inverting Debashis Paul asymptotic
				 * limit for the expected value of the eigenvalues.
				 */
//				float gamma = p_params.beta * (float)sPC/(float)p_nSimP;
				float gamma = (float)sPC/(float)p_nSimP;

				if (thres != 0.f)
				{
					rr += (i_mat.covEigVals[i] > thres * (1 + gamma) * sigma2) ? 1 : 0;
					i_mat.covEigVals[i] -= (i_mat.covEigVals[i] > thres * (1 + gamma) * sigma2) ? 
					                       std::min(i_mat.covEigVals[i], (1 + gamma) * sigma2) : 
					                       i_mat.covEigVals[i];
				}
				else
					i_mat.covEigVals[i] -= std::min(i_mat.covEigVals[i], (1 + gamma) * sigma2);
 #else
				i_mat.covEigVals[i] -= sigma2;
				i_mat.covEigVals[i]  = (i_mat.covEigVals[i] > thres * sigma2) ? 
				                        i_mat.covEigVals[i] : 0;
 #endif
#else // SPTWO || VARIANCE_THRESHOLD_1
				/* NOTE: The first step in SPTWO is a hard thresholding by beta*sigma2 */
				i_mat.covEigVals[i] = (i_mat.covEigVals[i] > sigma2) ? FLT_MAX : 0.f ;
#endif
				rank_variance  += i_mat.covEigVals[i];
			}

//			printf("\trr = %d\n",rr);

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (unsigned k = 0; k < r; ++k)
				i_mat.covEigVals[k] = 1.f / ( 1. + sigma2 / i_mat.covEigVals[k] );

			
			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Z' = X'*U
			productMatrix(i_mat.groupTranspose,
							  io_group[c],
							  i_mat.covEigVecs,
							  p_nSimP, r, sPC,
							  false, false);

			//! U * W
			float *eigVecs = i_mat.covEigVecs.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < sPC; ++i)
				*eigVecs++ *= i_mat.covEigVals[k];

			//! hX' = Z'*(U*W)'
			productMatrix(io_group[c],
							  i_mat.groupTranspose,
							  i_mat.covEigVecs,
							  p_nSimP, sPC, r,
							  false, true);

			//! Add baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group[c][k] += i_mat.baricenter[j];
		}
		else
		{
			//! rank = 0: set as baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group[c][k] = i_mat.baricenter[j];

			//! Avoid 0/0 in return statement
			total_variance = 1.f;
		}
	}

	// return percentage of captured variance
	return rank_variance / total_variance;
}

/**
 * @brief Implementation of computeBayesEstimateStep1_LR computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the full SVD of the data matrix using LAPACK.
 *
 * See computeBayesEstimateStep1_LR for information about the arguments.
 **/
float computeBayesEstimateStep1_LR_SVD_LAPACK(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;
#else
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch;
#endif
	const unsigned r    = p_params.rank;
	const unsigned mdim = std::min(sPC, p_nSimP);

	//! Variances
	float  rank_variance = 0.f;
	float total_variance = 0.f;

	for (unsigned c = 0; c < io_group.size(); c++)
	{
		//! Center group of patches
		centerData(io_group[c], i_mat.baricenter, p_nSimP, sPC);

		//! Compute SVD
		int info = matrixSVD(io_group[c], sPC, p_nSimP,
				i_mat.svd_S, i_mat.svd_VT, i_mat.svd_U,
				i_mat.svd_work, i_mat.svd_iwork);

		/* NOTE: What does matrixSVD return?
		 *
		 * matrixSVD assumes that matrices are stored by columns.
		 * If X is the data matrix stored by rows, matrixSVD computes the SVD of
		 * X' = V*S*U', and returns V and UT in column-major layout. If we consider
		 * them in row-major layout, we actually have U and VT. This explains
		 * the names of the variables.
		 *
		 * Therefore, in column-major layout, we have that:
		 *
		 * svd_VT - p_nSimP x mdim - columns are left  sing. vectors of X'
		 * svd_U  - mdim x sPC     - rows    are right sing. vectors of X'
		 *
		 */

		//! Substract sigma2 and compute variance captured by the r leading eigenvectors
		for (int i = 0; i < r; ++i)
		{
			// transform sing. val of noisy data matrix into eig. val of cov. matrix
			i_mat.svd_S[i] *= i_mat.svd_S[i]/(float)p_nSimP;
			i_mat.svd_S[i] -= std::min(i_mat.svd_S[i], sigma2);
			rank_variance  += i_mat.svd_S[i];
			total_variance += i_mat.svd_S[i];
		}

		for (int i = r; i < mdim; ++i)
		{
			i_mat.svd_S[i] *= i_mat.svd_S[i]/(float)p_nSimP;
			i_mat.svd_S[i] -= std::min(i_mat.svd_S[i], sigma2);
			total_variance += i_mat.svd_S[i];
		}

		/* NOTE: How do we perform the filtering?
		 *
		 * The filtering can be done by modifying the singular values only. Thus
		 * fX' = V*fS*U'.
		 */

		//! Apply Bayes' filter to singular values
		for (unsigned k = 0; k < r; ++k)
			i_mat.svd_S[k] = sqrtf( (float)p_nSimP * i_mat.svd_S[k] ) 
			               / ( 1. + sigma2 / i_mat.svd_S[k] );

		//! Multiply svd_S times svd_U or svd_VT depending which is smaller
		if (p_nSimP < sPC)
		{
			// Multiply first k cols of svd_VT (left sing. vectors of X')
			float *svdVT = i_mat.svd_VT.data();
			for (unsigned k = 0; k < r      ; ++k)
			for (unsigned i = 0; i < p_nSimP; ++i)
				*svdVT++ *= i_mat.svd_S[k];
		}
		else
		{
			// Multiply first k rows of svd_U (right sing. vectors of X')
			for (unsigned k = 0; k < r; ++k)
			{
				float *svdU = i_mat.svd_U.data() + k;
				for (unsigned i = 0; i < sPC; ++i, svdU += sPC)
					*svdU *= i_mat.svd_S[k];
			}
		}

		//! Multiply svd_VT*svd_S*svd_U = filtered(data)'
		productMatrix(io_group[c],
		              i_mat.svd_VT,
		              i_mat.svd_U,
		              p_nSimP, sPC, r,
		              false, false, true,
		              p_nSimP, mdim);

		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_group[c][k] += i_mat.baricenter[j];
	}

	// return percentage of captured variance
	return rank_variance / total_variance;
}

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
 *
 * @param io_group: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *		- groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
float computeBayesEstimateStep1_LR(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Check if all patches are equal
	/* In some images might have a dark noiseless frame (for instance the image
	 * might be the result of a transformation, and a black region correspond to
	 * the pixels that where before the frame before the registration). A group
	 * of patches with 0 variance create problems afterwards. Such groups have to
	 * be avoided.
	 * */

	/* TODO
	 * This fix only addresses the problem in which all patches are equal. A more
	 * subtle situation is when some patches are equal. These patches shouldn't be
	 * denoised, since it is extremely unlekely that such an event occurs.
	 * */

#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;
#else
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch;
#endif

	for (unsigned c = 0; c < io_group.size(); c++)
	for (unsigned j = 0; j < sPC; j++)
	{
		float v = io_group[c][j * p_nSimP];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_group[c][j * p_nSimP + i])
				goto not_equal;
	}

	//! All patches are equal ~ do nothing
	return 0.f;

not_equal:
	//! Not all patches are equal ~ denoise
	return
#if defined(USE_SVD_LAPACK)
		computeBayesEstimateStep1_LR_SVD_LAPACK(io_group, i_mat,
			io_nInverseFailed, p_params, p_nSimP);
#else
		computeBayesEstimateStep1_LR_EIG_LAPACK(io_group, i_mat,
			io_nInverseFailed, p_params, p_nSimP);
#endif

}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_groupNoisy: inputs all similar patches in the noisy image,
 *                         outputs their denoised estimates.
 * @param i_groupBasic: contains all similar patches in the basic image.
 * @param i_mat: contains :
 *    - groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the image;
 * @param p_params: see processStep2 for more explanations;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
float computeBayesEstimateStep2_FR(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float diagVal = p_params.beta * p_params.sigma * p_params.sigma;
#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime * p_imSize.channels;
#else
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_imSize.channels;
#endif

	//! Check if all patches are equal
	/* In some images might have a dark noiseless frame (for instance the image
	 * might be the result of a transformation, and a black region correspond to
	 * the pixels that where before the frame before the registration). A group
	 * of patches with 0 variance create problems afterwards. Such groups have to
	 * be avoided.
	 * */

	/* TODO
	 * This fix only addresses the problem in which all patches are equal. A more
	 * subtle situation is when some patches are equal. These patches shouldn't be
	 * denoised, since it is extremely unlekely that such an event occurs.
	 * */

	for (unsigned j = 0; j < sPC; j++)
	{
		float v = io_groupNoisy[j * p_nSimP];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_groupNoisy[j * p_nSimP + i])
				goto not_equal;
	}

	//! All patches are equal ~ do nothing
	return 0.f;

not_equal:
	//! Not all patches are equal ~ denoise

	//! Center 3D groups around their baricenter
	centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);
	centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC);

	//! Compute the covariance matrix of the set of similar patches
	covarianceMatrix(i_groupBasic, i_mat.covMat, p_nSimP, sPC);

	//! Bayes' Filtering
	for (unsigned k = 0; k < sPC; k++)
		i_mat.covMat[k * sPC + k] += diagVal;

	//! Compute the estimate
	if (inverseMatrix(i_mat.covMat, sPC) == EXIT_SUCCESS)
	{
		productMatrix(i_groupBasic, i_mat.covMat, io_groupNoisy, sPC, sPC, p_nSimP);
		for (unsigned k = 0; k < sPC * p_nSimP; k++)
			io_groupNoisy[k] -= diagVal * i_groupBasic[k];
	}
	else 
		io_nInverseFailed++;

	//! Add baricenter
	for (unsigned j = 0, k = 0; j < sPC; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
			io_groupNoisy[k] += i_mat.baricenter[j];

	return 1.f;
}

#if defined(FRAMES_DECOUPLED)
/**
 * @brief Implementation of computeBayesEstimateStep2_LR computing the
 * principal directions of the a priori covariance matrix. In this version
 * we consider that the different frames of the 3D patch are decoupled.
 * Thus instead of performing one Wiener filtering in IR^{sx sx st}
 * we do st independent Wiener filters in IR^{sx sx}.
 * This functions computes the eigenvectors/values of the data covariance
 * matrix using LAPACK.
 *
 * See computeBayesEstimateStep2_LR for information about the arguments.
 **/
float computeBayesEstimateStep2_LR_EIG_LAPACK(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	const bool p_flatPatch
){
	//! Parameters initialization
	const int sPx = p_params.sizePatch;
	const int sPt = p_params.sizePatchTime;
	const int chnls = p_imSize.channels;

	//! Alloc memory for the patch frames
	std::vector<float> frame_groupNoisy(sPx * sPx * chnls * p_nSimP);
	std::vector<float> frame_groupBasic(sPx * sPx * chnls * p_nSimP);

	float r_variance = 0.f;
	float total_variance = 1.f;

	for (int t = 0; t < sPt; ++t)
	{
		//! Adapt the data: store each patch frame in contiguous memory
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			unsigned kc = c * sPx * sPx * p_nSimP * sPt + 
				           t * sPx * sPx * p_nSimP;

			for (unsigned hy = 0; hy < sPx; hy++)
			for (unsigned hx = 0; hx < sPx; hx++)
			for (unsigned n = 0; n < p_nSimP; n++, k++, kc++)
			{
				frame_groupNoisy[k] = io_groupNoisy[kc];
				frame_groupBasic[k] =  i_groupBasic[kc];
			}
		}

		//! Wiener filtering of patch frame t starts here	
		const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
		const unsigned sPC = sPx * sPx * chnls;
		const unsigned r   = p_params.rank;

		//! Center 3D groups around their baricenter
		centerData(frame_groupBasic, i_mat.baricenter, p_nSimP, sPC);
		centerData(frame_groupNoisy, i_mat.baricenter, p_nSimP, sPC);

		if (r > 0)
		{
			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(frame_groupBasic, i_mat.covMat, p_nSimP, sPC);

			//! Compute total variance
			total_variance = 0.f;
			for (int i = 0; i < sPC; ++i)
				total_variance += i_mat.covMat[i*sPC + i];

			//! Compute leading eigenvectors
			int info = matrixEigs(i_mat.covMat, sPC, r, i_mat.covEigVals, i_mat.covEigVecs);

			//! Compute variance captured by the r leading eigenvectors
			for (int i = 0; i < r; ++i)
				r_variance += i_mat.covEigVals[i];

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (unsigned k = 0; k < r; ++k)
				i_mat.covEigVals[k] = 1.f / ( 1. + sigma2 / i_mat.covEigVals[k] );

			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Z' = X'*U
			productMatrix(i_mat.groupTranspose,
							  frame_groupNoisy,
							  i_mat.covEigVecs,
							  p_nSimP, r, sPC,
							  false, false);

			//! U * W
			float *eigVecs = i_mat.covEigVecs.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < sPC; ++i)
				*eigVecs++ *= i_mat.covEigVals[k];

			//! hX' = Z'*(U*W)'
			productMatrix(frame_groupNoisy,
							  i_mat.groupTranspose,
							  i_mat.covEigVecs,
							  p_nSimP, sPC, r,
							  false, true);

			//! Add baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					frame_groupNoisy[k] += i_mat.baricenter[j];
		}
		else
			//! r = 0: set all patches as baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					frame_groupNoisy[k] = i_mat.baricenter[j];

		//! Wiener filtering of patch frame t done

		//! Copy the filtered frame data back to 3D patches
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			unsigned kc = c * sPx * sPx * p_nSimP * sPt + 
				           t * sPx * sPx * p_nSimP;

			for (unsigned hy = 0; hy < sPx; hy++)
			for (unsigned hx = 0; hx < sPx; hx++)
			for (unsigned n = 0; n < p_nSimP; n++, k++, kc++)
				io_groupNoisy[kc] = frame_groupNoisy[k];
		}
	}

	// return percentage of captured variance
	return r_variance / total_variance;

}
#elif defined(CHANNELS_DECOUPLED)
/**
 * @brief Implementation of computeBayesEstimateStep2_LR computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the eigenvectors/values of the data covariance matrix using LAPACK.
 * In this version of the function, the color channels are treated
 * independently, as in the first step.
 *
 * See computeBayesEstimateStep2_LR for information about the arguments.
 **/
float computeBayesEstimateStep2_LR_EIG_LAPACK(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	const bool p_flatPatch
){
	//! Parameters initialization
	const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
	      float sigmab2 = p_params.beta * p_params.sigma_basic * p_params.sigma_basic;
	const float thres  = p_params.variThres;
#ifndef SPTWO
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
#else
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch;
#endif
	const unsigned r   = p_params.rank;
	const unsigned nch = p_imSize.channels;

	//! Center 3D groups around their baricenter
#ifdef NOISY_COVARIANCE
	//! Center noisy patches with their baricenter
	centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC*nch);

#else //NOISY_COVARIANCE
 #ifdef BARICENTER_BASIC
	//! Center basic patches with their baricenter
	centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC*nch);

	//! Remove basic's baricenter from Noisy patches
	for (unsigned j = 0, k = 0; j < sPC*nch; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
			io_groupNoisy[k] -= i_mat.baricenter[j];

 #else //BARICENTER_BASIC
	if (!p_flatPatch)
	{
		//! Center basic patches with their baricenter
		centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);

		//! Center noisy patches with their baricenter
		centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC);
	}
	else
	{
		//! Center basic patches with their baricenter
		centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);

		//! Remove basic's baricenter from Noisy patches
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] -= i_mat.baricenter[j];
	}

 #endif//BARICENTER_BASIC
#endif//NOISY_COVARIANCE

	float r_variance = 0.f;
	float total_variance = 1.f;

	if (r > 0)
	{
		for (unsigned c = 0; c < nch; c++)
		{

			std::vector<float> groupNoisy_c(io_groupNoisy.begin() + sPC*p_nSimP * c   ,
			                                io_groupNoisy.begin() + sPC*p_nSimP *(c+1));

			std::vector<float> groupBasic_c( i_groupBasic.begin() + sPC*p_nSimP * c   ,
			                                 i_groupBasic.begin() + sPC*p_nSimP *(c+1));

			//! Compute the covariance matrix of the set of similar patches
#ifndef NOISY_COVARIANCE
			covarianceMatrix(groupBasic_c, i_mat.covMat, p_nSimP, sPC);
#else
			covarianceMatrix(groupNoisy_c, i_mat.covMat, p_nSimP, sPC);
#endif

			//! Compute total variance
			total_variance = 0.f;
			for (int i = 0; i < sPC; ++i)
				total_variance += i_mat.covMat[i*sPC + i];

			//! Compute leading eigenvectors
			int info = matrixEigs(i_mat.covMat, sPC, r, i_mat.covEigVals, i_mat.covEigVecs);

#ifdef NOISY_COVARIANCE
			sigmab2 = sigma2;
#endif

			//! Convariance matrix is noisy, either because it has been computed
			//  from the noisy patches, or because we model the noise in the basic
			//  estimate
			if (sigmab2)
			{
				for (int i = 0; i < r; ++i)
				{
#ifdef CLIPPED_VARIANCE
					i_mat.covEigVals[i] -= std::min(i_mat.covEigVals[i], sigmab2);
#elif defined(PAUL_VARIANCE)
					/* variance estimate obtained by inverting Debashis Paul asymptotic
					 * limit for the expected value of the eigenvalues.
					 */
					float tmp, gamma = (float)sPC/(float)p_nSimP;

					if (i_mat.covEigVals[i] > sigmab2 * (tmp = (1 + sqrtf(gamma)))*tmp)
					{
						tmp = i_mat.covEigVals[i] - sigmab2 * (1 + gamma);
						i_mat.covEigVals[i] = tmp * 0.5 
						                    * (1. + sqrtf(std::max(0., 1. - 4.*gamma*sigmab2*sigmab2/tmp/tmp)));
					}
					else
						i_mat.covEigVals[i] = 0;

#elif defined(PAUL_SIMPLE_VARIANCE)
					float gamma = (float)sPC/(float)p_nSimP;
					i_mat.covEigVals[i] -= std::min(i_mat.covEigVals[i], (1 + gamma) * sigmab2);
#else
					i_mat.covEigVals[i] -= sigmab2;
#endif
				}
			}

			//! Compute variance captured by the r leading eigenvectors
			for (int i = 0; i < r; ++i)
				r_variance += i_mat.covEigVals[i];

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (unsigned k = 0; k < r; ++k)
#ifndef LI_ZHANG_DAI
				/* NOTE: thres is used to control a variance threshold proportional
				 * to the noise. If the eigenvalue is bellow the threshold, the 
				 * filter coefficient is set to zero. Otherwise it is set to 
				 * the clipped Wiener coefficient.
				 */
				i_mat.covEigVals[k] = (i_mat.covEigVals[k] > thres * sigma2) ?
				                      1.f / ( 1. + sigma2 / i_mat.covEigVals[k] ) :
				                      0.f;
#else
				i_mat.covEigVals[k] = (i_mat.covEigVals[k] > 4*sigma2)
				                    ? 0.5 * ( 1. + sqrt(1 - 4*sigma2 / i_mat.covEigVals[k])) 
				                    : 0;
#endif

			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Z' = X'*U
			productMatrix(i_mat.groupTranspose,
			              groupNoisy_c,
			              i_mat.covEigVecs,
			              p_nSimP, r, sPC,
			              false, false);

			//! U * W
			float *eigVecs = i_mat.covEigVecs.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < sPC; ++i)
				*eigVecs++ *= i_mat.covEigVals[k];

			//! hX' = Z'*(U*W)'
			productMatrix(groupNoisy_c,
			              i_mat.groupTranspose,
			              i_mat.covEigVecs,
			              p_nSimP, sPC, r,
			              false, true);

			//! Copy channel back into vector
			std::copy(groupNoisy_c.begin(), groupNoisy_c.end(),
			          io_groupNoisy.begin() + sPC*p_nSimP*c);
		}
	}
	else
	{
		//! r = 0: set all patches as baricenter
		for (unsigned j = 0, k = 0; j < sPC*nch; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] = 0;
	}

	//! Add baricenter
	for (unsigned j = 0, k = 0; j < sPC*nch; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] += i_mat.baricenter[j];


	// return percentage of captured variance
	return r_variance / total_variance;

}
#else
/**
 * @brief Implementation of computeBayesEstimateStep2_LR computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the eigenvectors/values of the data covariance matrix using LAPACK.
 *
 * See computeBayesEstimateStep2_LR for information about the arguments.
 **/
float computeBayesEstimateStep2_LR_EIG_LAPACK(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	const bool p_flatPatch
){
	//! Parameters initialization
	const float sigma2  = p_params.beta * p_params.sigma * p_params.sigma;
	      float sigmab2 = p_params.beta * p_params.sigma_basic * p_params.sigma_basic;
	const float thres   = p_params.variThres;
#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime * p_imSize.channels;
#else
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_imSize.channels;
#endif
	const unsigned r    = p_params.rank;

	//! Center 3D groups around their baricenter
#ifdef NOISY_COVARIANCE
	//! Center noisy patches with their baricenter
	centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC);

#else //NOISY_COVARIANCE
 #ifdef BARICENTER_BASIC
	//! Center basic patches with their baricenter
	centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);

	//! Remove basic's baricenter from Noisy patches
	for (unsigned j = 0, k = 0; j < sPC; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
			io_groupNoisy[k] -= i_mat.baricenter[j];

 #else //BARICENTER_BASIC
	if (!p_flatPatch)
	{
		//! Center basic patches with their baricenter
		centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);

		//! Center noisy patches with their baricenter
		centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC);
	}
	else
	{
		//! Center basic patches with their baricenter
		centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);

		//! Remove basic's baricenter from Noisy patches
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] -= i_mat.baricenter[j];
	}

 #endif//BARICENTER_BASIC
#endif//NOISY_COVARIANCE

	float r_variance = 0.f;
	float total_variance = 1.f;

	if (r > 0)
	{
		//! Compute the covariance matrix of the set of similar patches
#ifndef NOISY_COVARIANCE
		covarianceMatrix(i_groupBasic, i_mat.covMat, p_nSimP, sPC);
#else
		covarianceMatrix(io_groupNoisy, i_mat.covMat, p_nSimP, sPC);
#endif

		//! Compute total variance
		total_variance = 0.f;
		for (int i = 0; i < sPC; ++i)
			total_variance += i_mat.covMat[i*sPC + i];

		//! Compute leading eigenvectors
		int info = matrixEigs(i_mat.covMat, sPC, r, i_mat.covEigVals, i_mat.covEigVecs);

#ifdef NOISY_COVARIANCE
		sigmab2 = sigma2;
#endif

		//! Convariance matrix is noisy, either because it has been computed
		//  from the noisy patches, or because we model the noise in the basic
		//  estimate
		if (sigmab2)
		{
			for (int i = 0; i < r; ++i)
			{
#ifdef CLIPPED_VARIANCE
				i_mat.covEigVals[i] -= std::min(i_mat.covEigVals[i], sigmab2);
#elif defined(PAUL_VARIANCE)
				/* variance estimate obtained by inverting Debashis Paul asymptotic
				 * limit for the expected value of the eigenvalues.
				 */
				float tmp, gamma = (float)sPC/(float)p_nSimP;

				if (i_mat.covEigVals[i] > sigmab2 * (tmp = (1 + sqrtf(gamma)))*tmp)
				{
					tmp = i_mat.covEigVals[i] - sigmab2 * (1 + gamma);
					i_mat.covEigVals[i] = tmp * 0.5 
						                 * (1. + sqrtf(std::max(0., 1. - 4.*gamma*sigmab2*sigmab2/tmp/tmp)));
				}
				else
					i_mat.covEigVals[i] = 0;

#elif defined(PAUL_SIMPLE_VARIANCE)
				float gamma = (float)sPC/(float)p_nSimP;
				i_mat.covEigVals[i] -= std::min(i_mat.covEigVals[i], (1 + gamma) * sigmab2);
#else
				i_mat.covEigVals[i] -= sigmab2;
#endif
			}
		}


		//! Compute variance captured by the r leading eigenvectors
		for (int i = 0; i < r; ++i)
			r_variance += i_mat.covEigVals[i];

		//! Compute eigenvalues-based coefficients of Bayes' filter
		for (unsigned k = 0; k < r; ++k)
#ifndef LI_ZHANG_DAI
			/* NOTE: thres is used to control a variance threshold proportional
			 * to the noise. If the eigenvalue is bellow the threshold, the 
			 * filter coefficient is set to zero. Otherwise it is set to 
			 * the clipped Wiener coefficient.
			 */
			i_mat.covEigVals[k] = (i_mat.covEigVals[k] > thres * sigma2) ?
			                      1.f / ( 1. + sigma2 / i_mat.covEigVals[k] ) :
			                      0.f;
#else
			i_mat.covEigVals[k] = (i_mat.covEigVals[k] > 4*sigma2)
			                    ? 0.5 * ( 1. + sqrt(1 - 4*sigma2 / i_mat.covEigVals[k])) 
			                    : 0;
#endif

		/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
		 * row a patch. Thus, in column-major storage it corresponds to X^T, where
		 * each column of X contains a centered data point.
		 *
		 * We need to compute the noiseless estimage hX as 
		 * hX = U * W * U' * X
		 * where U is the matrix with the eigenvectors and W is a diagonal matrix
		 * with the filter coefficients.
		 *
		 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
		 * we compute 
		 * hX' = X' * U * (W * U')
		 */

		//! Z' = X'*U
		productMatrix(i_mat.groupTranspose,
		              io_groupNoisy,
		              i_mat.covEigVecs,
		              p_nSimP, r, sPC,
		              false, false);

		//! U * W
		float *eigVecs = i_mat.covEigVecs.data();
		for (unsigned k = 0; k < r  ; ++k)
		for (unsigned i = 0; i < sPC; ++i)
			*eigVecs++ *= i_mat.covEigVals[k];

		//! hX' = Z'*(U*W)'
		productMatrix(io_groupNoisy,
		              i_mat.groupTranspose,
		              i_mat.covEigVecs,
		              p_nSimP, sPC, r,
		              false, true);

		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] += i_mat.baricenter[j];
	}
	else
		//! r = 0: set all patches as baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] = i_mat.baricenter[j];


	// return percentage of captured variance
	return r_variance / total_variance;

}
#endif

/**
 * @brief Implementation of computeBayesEstimateStep2_LR computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the full SVD of the data matrix using LAPACK.
 *
 * See computeBayesEstimateStep2_LR for information about the arguments.
 **/
float computeBayesEstimateStep2_LR_SVD_LAPACK(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime * p_imSize.channels;
#else
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_imSize.channels;
#endif
	const unsigned r    = p_params.rank;
	const unsigned mdim = std::min(sPC, p_nSimP);

	//! Center 3D groups around their baricenter
	centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);
	centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC);

	//! Compute SVD
	int info = matrixSVD(i_groupBasic, sPC, p_nSimP,
			i_mat.svd_S, i_mat.svd_VT, i_mat.svd_U,
			i_mat.svd_work, i_mat.svd_iwork);

	//! Compute variance captured by the r leading eigenvectors
	float r_variance = 0.f;
	for (int i = 0; i < r; ++i)
	{
		i_mat.svd_S[i] *= i_mat.svd_S[i]/(float)p_nSimP;
		r_variance += i_mat.svd_S[i];
	}

	//! Compute total variance
	float total_variance = r_variance;
	for (int i = r; i < mdim; ++i)
	{
		i_mat.svd_S[i] *= i_mat.svd_S[i]/(float)p_nSimP;
		total_variance += i_mat.svd_S[i];
	}

	//! Compute eigenvalues-based coefficients of Bayes' filter
	for (unsigned k = 0; k < r; ++k)
		i_mat.svd_S[k] = 1.f / sqrtf( 1. + sigma2 / i_mat.svd_S[k] );

	//! Scale eigenvectors using the filter coefficients
	float *svdU = i_mat.svd_U.data();
	for (unsigned k = 0; k < r  ; ++k)
	for (unsigned i = 0; i < sPC; ++i)
		svdU[mdim*i + k] *= i_mat.svd_S[k];

	/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
	 * row a patch. Thus, in column-major storage it corresponds to X^T, where
	 * each column of X contains a centered data point.
	 *
	 * We need to compute the noiseless estimage hX as 
	 * hX = U * U' * X
	 * where U is the matrix with the normalized eigenvectors.
	 *
	 * Matrix U' is stored (column-major) in i_mat.svd_U. Since we have X^T
	 * we compute 
	 * hX' = X' * (U')' * U'
	 */

	//! Z' = X'*U
	productMatrix(i_mat.groupTranspose,
	              io_groupNoisy,
	              i_mat.svd_U,
	              p_nSimP, r, sPC,
	              false, true, true,
	              p_nSimP, mdim);

	//! hX' = Z'*U'
	productMatrix(io_groupNoisy,
	              i_mat.groupTranspose,
	              i_mat.svd_U,
	              p_nSimP, sPC, r,
	              false, false, true,
	              p_nSimP, mdim);

	//! Add baricenter
	for (unsigned j = 0, k = 0; j < sPC; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
			io_groupNoisy[k] += i_mat.baricenter[j];

	// return percentage of captured variance
	return r_variance / total_variance;
}

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
 *
 * @param io_groupNoisy: inputs all similar patches in the noisy image,
 *                         outputs their denoised estimates.
 * @param i_grougroup contains all similar patches in the basic image.
 * @param i_mat: contains :
 *    - groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the image;
 * @param p_params: see processStep2 for more explanations;
 * @param p_nSimP: number of similar patches.
 *
 * @return estimate of kept variance.
 **/
float computeBayesEstimateStep2_LR(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_size
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	const bool p_flatPatch
){
	//! Check if all patches are equal
	/* In some images might have a dark noiseless frame (for instance the image
	 * might be the result of a transformation, and a black region correspond to
	 * the pixels that where before the frame before the registration). A group
	 * of patches with 0 variance create problems afterwards. Such groups have to
	 * be avoided.
	 * */

	/* TODO
	 * This fix only addresses the problem in which all patches are equal. A more
	 * subtle situation is when some patches are equal. These patches shouldn't be
	 * denoised, since it is extremely unlekely that such an event occurs.
	 * */

#ifndef SPTWO
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime * p_size.channels;
#else
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_size.channels;
#endif

	for (unsigned j = 0; j < sPC; j++)
	{
		float v = io_groupNoisy[j * p_nSimP];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_groupNoisy[j * p_nSimP + i])
				goto not_equal;
	}

	//! All patches are equal ~ do nothing
	return 0.f;

not_equal:
	//! Not all patches are equal ~ denoise
	return
#if defined(USE_SVD_LAPACK)
		computeBayesEstimateStep2_LR_SVD_LAPACK(io_groupNoisy, i_groupBasic, i_mat,
			io_nInverseFailed, p_size, p_params, p_nSimP);
#else
		computeBayesEstimateStep2_LR_EIG_LAPACK(io_groupNoisy, i_groupBasic, i_mat,
			io_nInverseFailed, p_size, p_params, p_nSimP, p_flatPatch);
#endif

}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 **/
int computeAggregationStep1(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<std::vector<float> > const& i_group
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.sz.channels;
	const unsigned sPx   = p_params.sizePatch;
#ifndef SPTWO
	const unsigned sPt   = p_params.sizePatchTime;
#else
	const unsigned sPt   = 1;
#endif

	const unsigned w   = io_im.sz.width;
	const unsigned h   = io_im.sz.height;
	const unsigned wh  = io_im.sz.wh;
	const unsigned whc = io_im.sz.whc;

	int masked = 0;

#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;

		unsigned px0, py0, pt0, pc0;
		io_im.sz.coords(ind, px0, py0, pt0, pc0);

		// FIXME: this isn't well defined for 3D patches
		if (p_params.onlyFrame >= 0 && pt0 != p_params.onlyFrame)
			continue;

		for (unsigned pt = 0; pt < sPt; pt++)
		for (unsigned py = 0; py < sPx; py++)
		for (unsigned px = 0; px < sPx; px++)
#ifdef SPTWO // avoid aggregating simmetrized patches
		if ( px0 + px >= 0 && px0 + px < io_im.sz.width  && 
		     py0 + py >= 0 && py0 + py < io_im.sz.height && 
		     pt0 + pt >= 0 && pt0 + pt < io_im.sz.frames ) 
#endif
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * wh;
				io_im(ij + pt * whc + py * w + px) += 
					i_group[c][(pt * sPx*sPx + py * sPx + px) * p_nSimP + n];
			}
			io_weight(ind1 + pt * wh + py * w + px)++;
		}

		//! Use Paste Trick
		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (io_mask(ind1)) masked++;
		io_mask(ind1) = false;

		if (p_params.doPasteBoost)
		{
			if ((py >     2*sPx) && io_mask(ind1 - w)) masked++;
			if ((py < h - 2*sPx) && io_mask(ind1 + w)) masked++;
			if ((px >     2*sPx) && io_mask(ind1 - 1)) masked++;
			if ((px < w - 2*sPx) && io_mask(ind1 + 1)) masked++;

			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < h - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px < w - 2*sPx) io_mask(ind1 + 1) = false;
		}
	}
#else
	for (unsigned n = 0; n < 1; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0; pt < sPt; pt++)
		for (unsigned py = 0; py < sPx; py++)
		for (unsigned px = 0; px < sPx; px++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * wh;
				io_im(ij + pt * whc + py * w + px) += 
					i_group[c][(pt * sPx*sPx + py * sPx + px) * p_nSimP + n];
			}
			io_weight(ind1 + pt * wh + py * w + px)++;
		}
	}

	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;

		//! Use Paste Trick
		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (io_mask(ind1)) masked++;
		io_mask(ind1) = false;

		if (p_params.doPasteBoost)
		{
			if ((py >     2*sPx) && io_mask(ind1 - w)) masked++;
			if ((py < h - 2*sPx) && io_mask(ind1 + w)) masked++;
			if ((px >     2*sPx) && io_mask(ind1 - 1)) masked++;
			if ((px < w - 2*sPx) && io_mask(ind1 + 1)) masked++;

			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < h - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px < w - 2*sPx) io_mask(ind1 + 1) = false;
		}
	}
#endif

	return masked;
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 **/
int computeAggregationStep2(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<float> const& i_group
,	Video<float> &variance
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.sz.channels;
	const unsigned sPx   = p_params.sizePatch;
#ifndef SPTWO
	const unsigned sPt   = p_params.sizePatchTime;
#else
	const unsigned sPt   = 1;
#endif

	const unsigned w   = io_im.sz.width;
	const unsigned h   = io_im.sz.height;
	const unsigned wh  = io_im.sz.wh;
	const unsigned whc = io_im.sz.whc;

	//! Aggregate estimates
	int masked = 0;
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];

		unsigned px0, py0, pt0, pc0;
		io_im.sz.coords(ind, px0, py0, pt0, pc0);

		// FIXME: this isn't well defined for 3D patches
		if (p_params.onlyFrame >= 0 && pt0 != p_params.onlyFrame)
			continue;

		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned pt = 0; pt < sPt; pt++)
			for (unsigned py = 0; py < sPx; py++)
			for (unsigned px = 0; px < sPx; px++, k++)
#ifdef SPTWO // avoid aggregating simmetrized patches
			if ( px0 + px >= 0 && px0 + px < io_im.sz.width  && 
				  py0 + py >= 0 && py0 + py < io_im.sz.height && 
				  pt0 + pt >= 0 && pt0 + pt < io_im.sz.frames ) 
#endif
				io_im(ij + pt * whc + py * w + px) +=
					i_group[k * p_nSimP + n];
		}

		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0; pt < sPt; pt++)
		for (unsigned py = 0; py < sPx; py++)
		for (unsigned px = 0; px < sPx; px++)
#ifdef SPTWO // avoid aggregating simmetrized patches
		if ( px0 + px >= 0 && px0 + px < io_im.sz.width  && 
			  py0 + py >= 0 && py0 + py < io_im.sz.height && 
			  pt0 + pt >= 0 && pt0 + pt < io_im.sz.frames ) 
#endif
			io_weight(ind1 + pt * wh + py * w + px)++;

		//! Apply Paste Trick
		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (io_mask(ind1)) masked++;
		io_mask(ind1) = false;

		if (p_params.doPasteBoost)
		{
			if ((py >     2*sPx) && io_mask(ind1 - w)) masked++;
			if ((py < h - 2*sPx) && io_mask(ind1 + w)) masked++;
			if ((px >     2*sPx) && io_mask(ind1 - 1)) masked++;
			if ((px < w - 2*sPx) && io_mask(ind1 + 1)) masked++;

			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < h - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px < w - 2*sPx) io_mask(ind1 + 1) = false;

#ifdef DEBUG_SHOW_WEIGHT
			if (py >     2*sPx) variance(ind1 - w) = variance(ind1);
			if (py < h - 2*sPx) variance(ind1 + w) = variance(ind1);
			if (px >     2*sPx) variance(ind1 - 1) = variance(ind1);
			if (px < w - 2*sPx) variance(ind1 + 1) = variance(ind1);
#endif
		}
	}

	return masked;
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D
 * group. This version is for a test: in the original version, all patches
 * in the group are marked as processed, and cannot be origins of a patch
 * group. In this version we only mark as processed the patches of the 
 * group which are nearby frames to the group origin.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeTemporalAggregationStep1(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<std::vector<float> > const& i_group
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls  = io_im.sz.channels;
	const unsigned width  = io_im.sz.width;
	const unsigned height = io_im.sz.height;
	const unsigned sP     = p_params.sizePatch;

	//! Compute coordinates of group origin
	int FRAME_STEP = 2;
	unsigned ox, oy, ot, oc;
	io_im.sz.coords(i_index[0], ox, oy, ot, oc);
	
#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;
		for (unsigned p = 0; p < sP; p++)
		for (unsigned q = 0; q < sP; q++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * width * height;
				io_im(ij + p * width + q) += i_group[c][(p * sP + q) * p_nSimP + n];
			}
			io_weight(ind1 + p * width + q)++;
		}

		// compute coordinates of current pixel
		unsigned px, py, pt, pc;
		io_im.sz.coords(ind, px, py, pt, pc);

		if (abs((int)pt - (int)ot) < FRAME_STEP)
		{
			//! Use Paste Trick
			io_mask(ind1) = false;

			if (p_params.doPasteBoost)
			{
				io_mask(ind1 - width) = false;
				io_mask(ind1 + width) = false;
				io_mask(ind1 - 1    ) = false;
				io_mask(ind1 + 1    ) = false;
			}
		}
	}
#else
	for (unsigned n = 0; n < 1; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;
		for (unsigned p = 0; p < 1; p++)
		for (unsigned q = 0; q < 1; q++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * width * height;
				io_im(ij + p * width + q) += i_group[c][(p * sP + q) * p_nSimP + n];
			}
			io_weight(ind1 + p * width + q)++;
		}
	}

	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;

		// compute coordinates of current pixel
		unsigned px, py, pt, pc;
		io_im.sz.coords(ind, px, py, pt, pc);

		if (abs((int)pt - (int)ot) < FRAME_STEP)
		{
			//! Use Paste Trick
			io_mask(ind1) = false;

			if (p_params.doPasteBoost)
			{
				io_mask(ind1 - width) = false;
				io_mask(ind1 + width) = false;
				io_mask(ind1 - 1    ) = false;
				io_mask(ind1 + 1    ) = false;
			}
		}
	}
#endif
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D
 * group. This version is for a test: in the original version, all patches
 * in the group are marked as processed, and cannot be origins of a patch
 * group. In this version we only mark as processed the patches of the 
 * group which are nearby frames to the group origin.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeTemporalAggregationStep2(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<float> const& i_group
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.sz.channels;
	const unsigned width = io_im.sz.width;
	const unsigned wh    = width * io_im.sz.height;
	const unsigned sP    = p_params.sizePatch;

	//! Compute coordinates of group origin
	int FRAME_STEP = 2;
	unsigned ox, oy, ot, oc;
	io_im.sz.coords(i_index[0], ox, oy, ot, oc);
	
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind  = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned p = 0; p < sP; p++)
			for (unsigned q = 0; q < sP; q++, k++)
				io_im(ij + p * width + q) += i_group[k * p_nSimP + n];
		}

		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;
		for (unsigned p = 0; p < sP; p++)
		for (unsigned q = 0; q < sP; q++)
			io_weight(ind1 + p * width + q)++;

		// compute coordinates of current pixel
		unsigned px, py, pt, pc;
		io_im.sz.coords(ind, px, py, pt, pc);

		if (abs((int)pt - (int)ot) < FRAME_STEP)
		{
			//! Apply Paste Trick
			io_mask(ind1) = false;

			if (p_params.doPasteBoost)
			{
				io_mask(ind1 - width) = false;
				io_mask(ind1 + width) = false;
				io_mask(ind1 - 1    ) = false;
				io_mask(ind1 + 1    ) = false;
			}
		}
	}
}

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_im: image of reference, when the weight if null;
 * io_im: will contain the final image;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return : none.
 **/
void computeWeightedAggregation(
	Video<float> const& i_im
,	Video<float> &io_im
,	Video<float> const& i_weight
){
	for (unsigned f = 0; f < i_im.sz.frames  ; f++)
	for (unsigned c = 0; c < i_im.sz.channels; c++)
	for (unsigned y = 0; y < i_im.sz.height  ; y++)
	for (unsigned x = 0; x < i_im.sz.width   ; x++)
		if (i_weight(x,y,f) > 0.f) io_im(x,y,f,c) /= i_weight(x,y,f);
		else                       io_im(x,y,f,c) /= i_im(x,y,f,c);

}

}

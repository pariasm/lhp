/*
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work Copyright (c) 2014, Pablo Arias <parias@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef VIDEO_NL_BAYES_H_INCLUDED
#define VIDEO_NL_BAYES_H_INCLUDED

/* Compute the 2nd step covariance matrix from the noisy patches. In this way,
 * the basic estimate is used only in the computation of the patch distances.*/
//#define NOISY_COVARIANCE

/* Uses an adaptation of Li,Zhand,Dai fixed point iteration to estimate the
 * signal power in the empirical Wiener filter. It applies whenever the Gaussian
 * model is learnt from the noisy patches, but currently it is implemented only 
 * in the second step (it should always be used together with the NOISY_COVARIANCE
 * option defined). */
//#define LI_ZHANG_DAI

/* Use a 0-1 linear filter computed by thresholding the variances of the cov. matrix. */
//#define VARIANCE_THRESHOLD_1
//#define VARIANCE_THRESHOLD_2

/* Corrects the 'centering bug' discovered by Nicola. In the second step, basic
 * and noisy patches are centered using the basic baricenter. If left undefined,
 * each set of patches (noisy and basic) are centered using their own baricenter. */
//#define BARICENTER_BASIC

/* Select a variance estimation method.
 *
 * CLIPPED_VARIANCE
 * Avoids negative weights in the empirical Wiener filter. When the estimated
 * variance of a certain component is lower than the noise variance, the filter
 * coefficient is set to zero. This applies whenever the Gaussian model is
 * estimated from the noisy patches: in the first step, or in the second step
 * if NOISY_COVARIANCE option is defined.
 *
 * PAUL_VARIANCE
 * Uses the full inverse of the asymptotic limit of the expected eigenvalue
 * derived by Debashis Paul.
 *
 * PAUL_SIMPLE_VARIANCE
 * Uses a simplified inverse of the asymptotic limit of the expected eigenvalue
 * derived by Debashis Paul. 
 *
 * DEFAULT
 * Just substract sigma^2 from the diagonal.
 * */
#define CLIPPED_VARIANCE
//#define PAUL_VARIANCE
//#define PAUL_SIMPLE_VARIANCE

/* Decouple the frames of the 3D patches in the 2nd step. This implies that
 * each frame is considered independent of the others. Thus instead of
 * computing a single Gaussian model of dimensionality st*ch*sx*sx, we compute
 * st Gaussian models each of dimensionality ch*sx*sx. This greatly reduces the
 * computation time, but also the quality.*/
//#define FRAMES_DECOUPLED

/* Decouple the color channels of the 3D patches in the 2nd step. Each channel is
 * considered independently of the others as in the first step.*/
#define CHANNELS_DECOUPLED

/* Implementation of the SPTWO method of Buades, Lisani and Miladinović. Motion
 * compensated extended 3D patches are considered for distance computation, but
 * the Gaussian model is built from their 2D slices. */
//#define SPTWO

/* Flat area tricks : enable the original flat area trick. Otherwise default to
 * the one which uses the basic baricenter in the second step (and does nothing
 * in the first step)*/
//#define FAT_ORIGINAL


#include "../Utilities/LibVideoT.hpp"
#include "../Utilities/Utilities.h"

namespace VideoNLB
{

/**
 * @brief Structures of parameters dedicated to NL-Bayes process
 *
 * @param sigma: value of noise;
 * @param sizePatch: size of patches (sizePatch x sizePatch);
 * @param nSimilarPatches: minimum number of similar patches wanted;
 * @param sizeSearchWindow: size of the search window around the reference patch;
 * @param boundary: must be > sizeSearchWindow. Boundary kept around sub-images when the image is
 *        subdivided for parallelization;
 * @param offSet: step between two similar patches;
 * @param useHomogeneousArea: if true, use the homogeneous area trick;
 * @param gamma: threshold to detect homogeneous area;
 * @param beta: parameter used to estimate the covariance matrix;
 * @param tau: parameter used to determine similar patches;
 * @param isFirstStep: true if the first step of the algorithm is wanted;
 * @param doPasteBoost: if true, patches near denoised similar patches will not be used as reference
 *        patches;
 * @param verbose: if true, print informations.
 **/
struct nlbParams
{
	float sigma;
	float sigma_basic;         // model remanent noise in the basic estimate when estimating cov. matrix
	unsigned sizePatch;        // depends on sigma
	unsigned sizePatchTime;    // user given
	unsigned nSimilarPatches;  // depends on sigma, sizeSearchTimeRange (1 channel) or sizePatch (3 channels)
	unsigned sizeSearchWindow; // depends on nSimilarPatches
	unsigned sizeSearchTimeRangeFwd; // how many forward  frames in search cube
	unsigned sizeSearchTimeRangeBwd; // how many backward frames in search cube
	unsigned boundary;         // depends on sizeSearchWindow
	unsigned offSet;           // depends on sizePatch
	unsigned offSetTime;       // depends on sizePatchTime
	bool useHomogeneousArea;
	float gamma;
	float variThres;           // variance threshold 
	unsigned rank;             // rank of covariance matrix
	float beta;                // depends on sigma
	float tau;                 // depends on sizePatch
	bool isFirstStep;
	bool doPasteBoost;
	bool verbose;
	int onlyFrame;             // denoise only onlyFrame (-1 means process all frames)
};

/**
 * @brief Structure containing matrices used in the Bayesian estimation.
 *
 * @param groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 * @param baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 * @param covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 * @param covMatTmp: allocated memory. Used to process the Bayes estimate;
 * @param tmpMat: allocated memory. Used to process the Bayes estimate.
 **/
struct matWorkspace
{
	std::vector<float> groupTranspose;
	std::vector<float> baricenter;
	std::vector<float> covMat;
	std::vector<float> covMatTmp;
	std::vector<float> tmpMat;

	// the following have been added for computing a low rank
	// approximation of the covariance matrix C

	// used if low-rank approximation is done via eigendecomp. of C
	std::vector<float> covEigVecs;
	std::vector<float> covEigVals;

	// used if low-rank approximation of C is done via SVD of data matrix X
	std::vector<float> svd_U;      // left  sing. vecs of X ~ eigen vecs of C
	std::vector<float> svd_V;      // right sing. vecs of X
	std::vector<float> svd_VT;     // right sing. vecs of X (transposed for LAPACKE)
	std::vector<float> svd_S;      // sing. values of X ~ sqrt of eigen vecs of C
	std::vector<float> svd_work;
	std::vector<int  > svd_iwork;

	// double matrix to compute svd using idd
	std::vector<double> svd_dU;    // left sing. vecs of X ~ eigen vecs of C
	std::vector<double> svd_dV;    // right sing. vecs of X
	std::vector<double> svd_dS;    // sing. values of X ~ sqrt of eigen vecs of C
	std::vector<double> svd_ddata; // 
	std::vector<double> svd_dwork;
};

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
,	const bool p_useArea
,	const bool p_verbose
,	const unsigned timeSearchRangeFwd = 0
,	const unsigned timeSearchRangeBwd = 0
,	const unsigned sizePatchTime = 1
,	const unsigned rank = 4
);


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
void setSizeSearchWindow(
	nlbParams& prms
,	unsigned sizeSearchWindow
);


/**
 * @brief Sets size of the patch. It sets the pixel offset as half the patch
 * size (this is BM3D speed-up).
 *
 * @param prms      : nlbParams for first or second step of the algorithm;
 * @param sizePatch : size of the patch;
 *
 * @return none.
 **/
void setSizePatch(
	nlbParams& prms
,	const VideoSize &p_size
,	unsigned sizePatch
);


/**
 * @brief Sets number of similar patches, ensuring that the number of similar
 * patches is not larger that the total number of available patches.
 *
 * @param prms            : nlbParams for first or second step of the algorithm;
 * @param nSimilarPatches : number of similar patches;
 *
 * @return none.
 **/
void setNSimilarPatches(
	nlbParams& prms
, unsigned nSimilarPatches
);


/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_params : nlbParams for first or second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	const nlbParams &i_params
);

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy video;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> &i_imNoisy
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,	const bool p_verbose
,  Video<float> &i_imClean
);

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy video;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_params1 : parameters for first step
 * @param p_params1 : parameters for second step
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> &i_imNoisy
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_params1
,	const nlbParams p_params2
,  Video<float> &i_imClean
);

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy video;
 * @param i_fflow  : forward optical flow;
 * @param i_bflow  : backward optical flow;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_params1 : parameters for first step
 * @param p_params1 : parameters for second step
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> & i_imNoisy
,  Video<float> const& i_fflow
,  Video<float> const& i_bflow
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_params1
,	const nlbParams p_params2
,  Video<float> &i_imClean
);

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy video;
 * @param i_fflow: optical flow
 * @param i_bflow: optical flow
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: parameters of the method, contains:
 *			- sigma: standard deviation of the noise;
 *			- sizePatch: size of patches (sizePatch x sizePatch);
 *			- nSimilarPatches: number of similar patches;
 *			- sizeSearchWindow: size of the neighbourhood searching window;
 *			- useHomogeneousArea: if true, the trick of using homogeneous area will be used;
 *			- gamma: parameter used to determine if we are in an homogeneous area;
 *			- maxAvoid: parameter used to stop the paste trick;
 *			- beta: parameter used during the estimate of the denoised patch;
 *			- coefBaricenter: parameter to determine if the covariance matrix inversion is correct;
 *			- isFirstStep: true if it's the first step of the algorithm which is needed;
 *			- verbose: if true, print some informations, do nothing otherwise.
 *
 * @return Percentage of processed groups over number of pixels.
 **/
unsigned processNlBayes(
	Video<float> const& i_imNoisy
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	Video<float> &io_imBasic
,	Video<float> &o_imFinal
,	nlbParams const& p_params
,	VideoUtils::CropPosition p_crop = VideoUtils::CropPosition()
,  Video<float> const &i_imClean = Video<float>()
);

/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param i_im: contains the noisy video on which distances are processed;
 * @param o_group: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep1(
	Video<float> const& i_im
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	std::vector<std::vector<float> > &o_group
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const nlbParams &p_params
,  Video<float> const &i_imClean = Video<float>()
);

/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy video;
 * @param i_imBasic: contains the basic estimation;
 * @param o_groupNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_groupBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video<float> const& i_imNoisy
,	Video<float> const& i_imBasic
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	std::vector<float> &o_groupNoisy
,	std::vector<float> &o_groupBasic
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const nlbParams &p_params
,  Video<float> const &i_imClean = Video<float>()
);

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group: contains for each channels values of similar patches. If an homogeneous area
 *			is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the video.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep1(
	std::vector<std::vector<float> > &io_group
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
);

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
	std::vector<float> & io_groupNoisy
,	std::vector<float> const &i_groupBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
);

/**
 * @brief Compute the Bayes estimation.
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
float computeBayesEstimateStep1_FR(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
);

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
);

/**
 * @brief Compute the Bayes estimation.
 *
 * @param i_groupNoisy: contains all similar patches in the noisy video;
 * @param io_groupBasic: contains all similar patches in the basic video. Will contain estimates
 *			for all similar patches;
 * @param i_mat: contains :
 *		- groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the video;
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
);

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
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
float computeBayesEstimateStep2_LR(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,  const bool p_flatPatch = false
);

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the video with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
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
,	const nlbParams& p_params
,	const unsigned p_nSimP
);

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the video with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 *
 **/
int computeAggregationStep2(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<float> const& i_group
,	Video<float> &variance
,	std::vector<unsigned> const& i_index
,	const nlbParams& p_params
,	const unsigned p_nSimP
);

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
);

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
);

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_imReference: video of reference, when the weight if null;
 * io_imResult: will contain the final video;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return none.
 **/
void computeWeightedAggregation(
	Video<float> const& i_im
,	Video<float> &io_im
,	Video<float> const& i_weight
);

} // namespace

#endif // VIDEO_NL_BAYES_H_INCLUDED

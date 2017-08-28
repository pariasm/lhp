/**
 * \file imdiff.c
 * \brief Image difference calculator program
 * \author Pascal Getreuer <getreuer@cmla.ens-cachan.fr>
 *
 * This file implements the imdiff program, a command line tool for comparing
 * two images with various image quality metrics.  The possible metrics are
 * 
 *    - Max absolute difference, \f$ \max_n |A_n - B_n| \f$
 *    - Mean squared error, \f$ \tfrac{1}{N} \sum |A_n - B_n|^2 \f$
 *    - Root mean squared error, \f$ \sqrt{\mathrm{MSE}} \f$
 *    - Peak signal-to-noise ratio, \f$ -10 \log_{10}(\mathrm{MSE}/255^2) \f$
 * 
 * where N is the number of pixels and subscript n denotes the nth pixel.
 *
 * The program can also create a difference image, computed as
 * \f[ D_n = 255/20 (A_n - B_n) + 255/2 \f]
 * where values outside of the range [0,255] are saturated.
 *
 * \note Image alpha channels are ignored.  Also beware that although the
 * program can read 16-bit PNG images (provided USE_LIBPNG is enabled), the
 * image data is quantized internally to 8 bits.
 *
 * Copyright (c) 2010-2013, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version, or the terms of the
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <ctype.h>
#include <math.h>
#include <string.h>
#include "imageio.h"

/** \brief Display metrics for intensities in the range [0,DISPLAY_SCALING]. */
#define DISPLAY_SCALING     255

/** \brief enum of possible metrics */
typedef enum {DEFAULT_METRICS, MAX_METRIC, MSE_METRIC, RMSE_METRIC,
    PSNR_METRIC} metric_t;

/** \brief struct of program parameters. */
typedef struct
{
    /** \brief Input file A (clean). */
    char *file_a;
    /** \brief Input file B (distorted). */
    char *file_b;
    /** \brief Quality for saving JPEG images (0 to 100). */
    int jpeg_quality;
    /** \brief metric */
    metric_t metric;
    /** \brief Compute metric separately for each channel. */
    int separate_channels;
    /** \brief Ignore boundary effects by shaving a margin of size pad. */
    int pad;
        
    /** \brief Difference file. */
    char *difference_file;
    /** \brief Parameter D for creating the difference image. */
    float D;
} programparams;
    

void make_difference_image(float *A, const float *B,
    int width, int height, int num_channels, float D);
void basic_metrics(float *max_diff, float *mse, const float *A, const float *B,
    int width, int height, int num_channels, int pad);
int parse_params(programparams *param, int argc, char *argv[]);

/** \brief Print program usage help message. */
void print_help_message()
{
    puts("Image difference calculator, P. Getreuer 2010-2011, 2013\n");
    puts("Syntax: imdiff [options] <exact file> <distorted file>\n"
        "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:");
    puts("   -m <metric>  metric to use for comparison, choices are");
    puts("        max     Max absolute difference, max_n |A_n - B_n|");
    puts("        mse     Mean squared error, 1/N sum |A_n - B_n|^2");
    puts("        rmse    Root mean squared error, (MSE)^1/2");
    puts("        psnr    Peak signal-to-noise ratio, -10 log10(MSE/255^2)");
    puts("   -s           Compute metric separately for each channel");
    puts("   -p <pad>     Remove a margin of <pad> pixels before comparison");
    puts("   -D <number>  D parameter for difference image\n");
#ifdef USE_LIBJPEG
    puts("   -q <number>  Quality for saving JPEG images (0 to 100)\n");
#endif
    puts("Alternatively, a difference image is generated by the syntax\n"
    "   imdiff [-D <number>] <exact file> <distorted file> <output file>\n");
    puts("The difference image is computed as\n"
           "   D_n = 255/D (A_n - B_n) + 255/2.\n"
           "Values outside of the range [0,255] are saturated.\n");
    puts("Example:\n"
#ifdef USE_LIBPNG
        "   imdiff -mpsnr frog-exact.png frog-4x.png");
#else
        "   imdiff -mpsnr frog-exact.bmp frog-4x.bmp");
#endif
}

int main(int argc, char *argv[])
{
    struct
    {
        float *data;
        int width;
        int height;
    } A = {NULL, 0, 0}, B = {NULL, 0, 0};
    programparams param;
    float max_diff, max_diff_c[3], mse, mse_c[3];
    int channel, status = 1;
           
    if(!parse_params(&param, argc, argv))
        return 0;
    
    /* Read the exact image. */
    if(!(A.data = (float *)read_image(&A.width, &A.height, param.file_a,
        IMAGEIO_FLOAT | IMAGEIO_RGB | IMAGEIO_PLANAR)))
        goto fail;
    
    /* Read the distorted image. */
    if(!(B.data = (float *)read_image(&B.width, &B.height, param.file_b,
        IMAGEIO_FLOAT | IMAGEIO_RGB | IMAGEIO_PLANAR)))
        goto fail;
    
    if(A.width != B.width || A.height != B.height)
    {
        fprintf(stderr, "Image sizes don't match, %dx%d vs. %dx%d.\n",
            A.width, A.height, B.width, B.height);
        goto fail;
    }
    else if(A.width <= 2 * param.pad || A.height <= 2 * param.pad)
    {
        fprintf(stderr,
            "Removal of %d-pixel padding removes entire %dx%d image.\n",
            param.pad, A.width, A.height);
        goto fail;
    }
    
    if(param.difference_file)
    {
        make_difference_image(A.data, B.data, A.width, A.height, 3, param.D);
        
        if(!(write_image(A.data, A.width, A.height, param.difference_file,
            IMAGEIO_FLOAT | IMAGEIO_RGB | IMAGEIO_PLANAR, param.jpeg_quality)))
            goto fail;
    }
    else
    {
        max_diff = 0.0f;
        mse = 0.0f;
        
        for(channel = 0; channel < 3; channel++)
        {
            basic_metrics(&max_diff_c[channel], &mse_c[channel],
                          A.data + channel*A.width*A.height,
                          B.data + channel*B.width*B.height,
                          A.width, A.height, 1, param.pad);
           
            if(max_diff_c[channel] > max_diff)
                max_diff = max_diff_c[channel];
            
            mse += mse_c[channel];
        }
        
        mse /= 3;
        
        switch(param.metric)
        {
            case DEFAULT_METRICS:
                if(!param.separate_channels)
                {
                    printf("Maximum absolute difference:  %g\n",
                           DISPLAY_SCALING * max_diff);
                    printf("Root mean squared error:      %.4f\n", 
                           DISPLAY_SCALING * sqrt(mse));
                    printf("Peak signal-to-noise ratio:   %.4f\n",
                           -10 * log10(mse));
                }
                else
                {
                    printf("Maximum absolute difference:  %g %g %g\n",
                           DISPLAY_SCALING * max_diff_c[0],
                           DISPLAY_SCALING * max_diff_c[1],
                           DISPLAY_SCALING * max_diff_c[2]);
                    printf("Root mean squared error:      %.4f %.4f %.4f\n",
                           DISPLAY_SCALING * sqrt(mse_c[0]),
                           DISPLAY_SCALING * sqrt(mse_c[1]),
                           DISPLAY_SCALING * sqrt(mse_c[2]));
                    printf("Peak signal-to-noise ratio:   %.4f %.4f %.4f\n",
                           -10 * log10(mse_c[0]),
                           -10 * log10(mse_c[1]),
                           -10 * log10(mse_c[2]));
                }
                break;
            case MAX_METRIC:
                if(!param.separate_channels)
                    printf("%g\n", DISPLAY_SCALING*max_diff);
                else
                    printf("%g %g %g\n",
                           DISPLAY_SCALING * max_diff_c[0],
                           DISPLAY_SCALING * max_diff_c[1],
                           DISPLAY_SCALING * max_diff_c[2]);
                break;
            case MSE_METRIC:
                if(!param.separate_channels)
                    printf("%.4f\n", DISPLAY_SCALING * DISPLAY_SCALING * mse);
                else
                    printf("%.4f %.4f %.4f\n",
                           DISPLAY_SCALING * DISPLAY_SCALING * mse_c[0],
                           DISPLAY_SCALING * DISPLAY_SCALING * mse_c[1],
                           DISPLAY_SCALING * DISPLAY_SCALING * mse_c[2]);
                break;
            case RMSE_METRIC:
                if(!param.separate_channels)
                    printf("%.4f\n", DISPLAY_SCALING * sqrt(mse));
                else
                    printf("%.4f %.4f %.4f\n",
                           DISPLAY_SCALING * sqrt(mse_c[0]),
                           DISPLAY_SCALING * sqrt(mse_c[1]),
                           DISPLAY_SCALING * sqrt(mse_c[2]));
                break;
            case PSNR_METRIC:
                if(!param.separate_channels)
                    printf("%.4f\n", -10 * log10(mse));
                else
                    printf("%.4f %.4f %.4f\n",
                          -10 * log10(mse_c[0]),
                          -10 * log10(mse_c[1]),
                          -10 * log10(mse_c[2]));
                break;
        }
    }
    
    status = 0;
fail:
    if(B.data)
        free(B.data);
    if(A.data)
        free(A.data);
    return status;
}

/** \brief Make a difference image, diff = (A - B) / D + 0.5. */
void make_difference_image(float *A, const float *B,
    int width, int height, int num_channels, float D)
{
    const int num_el = num_channels * width * height;
    int n;
    
    D /= 255;
    
    for(n = 0; n < num_el; n++)
        A[n] = (A[n] - B[n]) / D + 0.5f;
}

/** \brief Compute the maximum absolute difference and the MSE. */
void basic_metrics(float *max_diff, float *mse, const float *A, const float *B,
    int width, int height, int num_channels, int pad)
{
    float diff, cur_max = 0;
    double accum_mse = 0;
    int x, y, channel, n;
    
    for(channel = 0; channel < num_channels; channel++)
        for(y = pad; y < height - pad; y++)
            for(x = pad; x < width - pad; x++)
            {
                n = x + width * (y + height * channel);
                diff = (float)fabs(A[n] - B[n]);
                    
                if(cur_max < diff)
                    cur_max = diff;
                
                accum_mse += diff*diff;
            }
    
    *max_diff = cur_max;
    *mse = (float)(
        accum_mse / (num_channels * (width - 2 * pad) * (height - 2 * pad)));
}

/** \brief Parse parameters from command line arguments. */
int parse_params(programparams *param, int argc, char *argv[])
{
    char *option_string;
    char option_char;
    int i;

    if(argc < 2)
    {
        print_help_message();
        return 0;
    }
    
    /* Set parameter defaults. */
    param->file_a = NULL;
    param->file_b = NULL;
    param->metric = DEFAULT_METRICS;
    param->separate_channels = 0;
    
    param->pad = 0;
    param->difference_file = NULL;
    param->jpeg_quality = 95;
    param->D = 20;
        
    for(i = 1; i < argc;)
    {
        if(argv[i] && argv[i][0] == '-')
        {
            if((option_char = argv[i][1]) == 0)
            {
                fprintf(stderr, "Invalid parameter format.\n");
                return 0;
            }

            if(argv[i][2])
                option_string = &argv[i][2];
            else if(++i < argc)
                option_string = argv[i];
            else
            {
                fprintf(stderr, "Invalid parameter format.\n");
                return 0;
            }
            
            switch(option_char)
            {
            case 'p':
                param->pad = atoi(option_string);
                
                if(param->pad < 0)
                {
                    fprintf(stderr, "pad must be nonnegative.\n");
                    return 0;
                }
                break;
            case 's':
                param->separate_channels = 1;
                i--;
                break;
            case 'D':
                param->D = (float)atof(option_string);

                if(param->D <= 0)
                {
                    fprintf(stderr, "D must be positive.\n");
                    return 0;
                }
                break;
            case 'm':
                if(!strcmp(option_string, "max"))
                    param->metric = MAX_METRIC;
                else if(!strcmp(option_string, "mse"))
                    param->metric = MSE_METRIC;
                else if(!strcmp(option_string, "rmse"))
                    param->metric = RMSE_METRIC;
                else if(!strcmp(option_string, "psnr"))
                    param->metric = PSNR_METRIC;
                else
                    fprintf(stderr, "Unknown metric.\n");
                break;
                
#ifdef USE_LIBJPEG
            case 'q':
                param->jpeg_quality = atoi(option_string);

                if(param->jpeg_quality <= 0 || param->jpeg_quality > 100)
                {
                    fprintf(stderr,
                            "JPEG quality must be between 0 and 100.\n");
                    return 0;
                }
                break;
#endif
            case '-':
                print_help_message();
                return 0;
            default:
                if(isprint(option_char))
                    fprintf(stderr, "Unknown option \"-%c\".\n", option_char);
                else
                    fprintf(stderr, "Unknown option.\n");

                return 0;
            }

            i++;
        }
        else
        {
            if(!param->file_a)
                param->file_a = argv[i];
            else if(!param->file_b)
                param->file_b = argv[i];
            else
                param->difference_file = argv[i];

            i++;
        }
    }
    
    if(!param->file_a || !param->file_b)
    {
        print_help_message();
        return 0;
    }
    
    return 1;
}

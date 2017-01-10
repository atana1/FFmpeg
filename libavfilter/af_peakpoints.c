/*
 * Copyright (c) 2016 Atana
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "libavcodec/avcodec.h"
#include "libavcodec/avfft.h"
#include "libavformat/avformat.h"
#include "libswscale/swscale.h"
#include "avfilter.h"
#include "audio.h"
#include "libavutil/opt.h"

#define SIZECHECK 4096
/* Upper and Lower Limit of sound frequency in audible range */
#define UPPER_LIMIT 300
#define LOWER_LIMIT 20
/* Total number of points in a single fingerprint, each point comes from a bin*/
#define KEYPOINTCOUNT 5


/* structure to store key points and time of its appearence */
typedef struct {
    size_t time;
    // array to store key points
    int *keyPoints;
} FingerPrint;


/* Structure to contain peak points context */
typedef struct {
    const AVClass *class;
    double *data;
    int nsamples;
    int index;
    size_t time;
    int isOnce;
    int buffFlag;
    int bin_size;
    FingerPrint *fprints;
    //double *peaks;
    int size; // number of peaks
    int windowSize;
    //char *inputFile;
} PeakPointsContext;


/* returns maximum value from an array conditioned on start and end index */
static double getMax(double *res_arr, int startIndex, int endIndex) {
    int i;
    double max = res_arr[startIndex];
    for (i = startIndex; i <= endIndex; i++) {
	    if (res_arr[i] > max) {
	        max = res_arr[i];
	    }
    }
    return max;
}

/* Stores peak frequency for each window(of chunkSize) in peaks array */
static void getPeakPointInChunk(int chunkSize, double *res_arr, int size, double *peaks) {
    int i = 0, peakIndex = 0;
    int startIndex = 0;
    double max;
    // get a chunk and find max value in it
    while (i < size) {
	    if (i % chunkSize-1 == 0) {
            max = getMax(res_arr, startIndex, i);
            peaks[peakIndex++] = max;
            startIndex = startIndex + chunkSize;
       }
       i += 1;
    }
}

/* Find out the range */
static int getIndex(int freq, int *range) {
    int i = 0;
    void *avc;

    // null check
    if (!range) {
        av_log(avc, AV_LOG_ERROR, "range is NULL\n");
    }

    while (range[i] < freq) {
        i++;
    }

    return i;
}

/* absolute value of complex number */
static double getAbs(FFTComplex cn) {
    double res;

    res = pow(cn.re, 2) + pow(cn.im, 2);
    res = sqrt(res);
    return res;
}

/* Finger printing */
static void getFingerPrint(FFTComplex *tab, int *bins, int bin_size, size_t time,
        FingerPrint *fp) {
    void *avc;
    int freq, index;
    double mag;
    int *range, *highscores;

    // allocate memmory for fp if not allocated
    if (!fp) {
        fp = av_malloc(sizeof(FingerPrint));
        fp->keyPoints = av_malloc_array(KEYPOINTCOUNT, sizeof(int));
    }

    range = av_malloc_array(KEYPOINTCOUNT, sizeof(int));

    // memmory check
    if (!range) {
        av_log(avc, AV_LOG_ERROR, "Can't allocate memmory for range array\n");
    }

    memcpy(range, bins, bin_size * sizeof(int));

    highscores = av_malloc_array(KEYPOINTCOUNT, sizeof(int));

    // memmory check
    if (!highscores) {
        av_log(avc, AV_LOG_ERROR, "Can't allocate memmory for highscore array\n");
    }

    for (freq = LOWER_LIMIT; freq < UPPER_LIMIT-1; freq++) {
        // Get magnitude, in decibel
        mag = log(getAbs(tab[freq]) + 1);

        // Find out the range in which it lies
        index = getIndex(freq, range);

        // Save the highest magnitude and corresponding frequency
        if (mag > highscores[index]) {
            highscores[index] = mag;
            fp->keyPoints[index] = freq;
        }
    }

    // key point from each bins are collected; assign time
    fp->time = time;
}

/* Get peaks points from windowed frequency domain data*/
static int getPeakPoints(PeakPointsContext *ppc) {
    int i, m, k, size, chunkSize, pSize, chunkSampleSize, resSize, fft_size;
    int j, inverse;
    double *fft_res;
    void *avc;
    FingerPrint fp;
    //RDFTContext *rdftC;
    FFTContext *fftc;
    FFTComplex *tab;
    //FFTSample *data;

    // frequency bins
    int bins[] = {40, 80, 120, 180, UPPER_LIMIT-1};
    int bin_size = 5;
    ppc->bin_size = bin_size;

    size = ppc->index;
    m = log2(ppc->windowSize);
    fft_size = 1 << m;
    chunkSize = ppc->windowSize;
    chunkSampleSize = size/chunkSize;
    resSize = chunkSize * chunkSampleSize;

    //fft_res = av_malloc_array(resSize, sizeof(double));

    /*if (!fft_res) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memmory for storing fft data\n");
        return 0;
    }*/

    ppc->fprints = av_malloc_array(resSize, sizeof(FingerPrint));

    // memmory check
    if (!ppc->fprints) {
        av_log(avc, AV_LOG_ERROR, "Can't allocate memmory for fingerprints\n");
    }

    // fft_size or chunkSize
    tab = av_malloc(chunkSize * sizeof(FFTComplex));

    // memmory check
    if (!tab) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memmory for complex data\n");
    }

    //rdftC = av_rdft_init(m, DFT_R2C);
    inverse = 0;
    fftc = av_fft_init(m, inverse);

    /*data = av_malloc_array(chunkSize, sizeof(FFTSample));

    if (!data) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memmory for chunk fft data\n");
        return 0;
    }*/

    // FFT transform for windowed time domain data
    // window is of size chunkSize
    k = 0;
    j = 0;
    while (k < resSize) {
        //copy data; put imaginary part as 0
        for (i = 0; i < chunkSize; i++) {
            //data[i] = ppc->data[i+k];
            tab[i].re = ppc->data[i+k];
            tab[i].im = 0.0;
        }
        //calculate FFT
        //av_rdft_calc(rdftC, data);
        av_fft_permute(fftc, tab);
        av_fft_calc(fftc, tab);

        /*for (i = 0; i < chunkSize; i++) {
            fft_res[i+k] = data[i];
        }*/

        // get finger print
        getFingerPrint(tab, bins, bin_size, ppc->time, &fp);
        ppc->fprints[j] = fp;

        j = j + 1;
        k = k + chunkSize;
    }

    //av_rdft_end(rdftC);
    av_fft_end(fftc);
    pSize = resSize/chunkSize;
    ppc->size = pSize;
    //ppc->peaks = av_malloc_array(pSize, sizeof(double));

    /*if (!ppc->peaks) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memory for peak storage\n");
        return 0;
    }*/

    //getPeakPointInChunk(chunkSize, fft_res, resSize, ppc->peaks);
    //av_freep(&data);
    av_freep(&tab);
    //av_freep(&fft_res);
    return 1;
}


#define OFFSET(x) offsetof(PeakPointsContext, x)

static const AVOption peakpoints_options[] = {
    { "wsize",  "set window size", OFFSET(windowSize),  AV_OPT_TYPE_INT,    {.i64=16},    0, INT_MAX},
    { NULL },
};

AVFILTER_DEFINE_CLASS(peakpoints);

static av_cold int init(AVFilterContext *ctx)
{
    PeakPointsContext *p = ctx->priv;

    if (p->windowSize < 16 || p->windowSize > SIZECHECK) {
        av_log(ctx, AV_LOG_ERROR, "window size must be in range 16 to %d\n", SIZECHECK);
        return AVERROR(EINVAL);
    }

    p->index = 0;
    p->time = 0;
    p->size = 0;
    p->isOnce = 1;
    p->data = av_malloc_array(SIZECHECK, sizeof(double));

    if (!p->data) {
        av_log(ctx, AV_LOG_ERROR, "Cann't allocate memmory for audio data\n");
        return AVERROR(EINVAL);
    }

    return 0;
}

static void ppointsStats(AVFilterContext *ctx, PeakPointsContext *p) {
    int i, j, ret;
    ret = getPeakPoints(p);

    if (ret && p->size) {
    	// print peaks
        if (p->isOnce) {
            av_log(ctx, AV_LOG_INFO, "######## Finger prints are ########\n");
            p->isOnce = 0;
        }

        for (i = 0; i < p->size; i++) {
            for (j = 0; j < p->bin_size; j++) {
                av_log(ctx, AV_LOG_INFO, "%d ", p->fprints[i].keyPoints[j]);
            }
            av_log(ctx, AV_LOG_INFO, "at timestep: %d\n", p->fprints[i].time);
        }

        // free peaks and set size to zero
        //av_freep(&p->peaks);
        av_freep(&p->fprints); // check if memory to be freed hierarchial fashion
        p->size = 0;
    } else if (p->size || !ret) {
        av_log(ctx, AV_LOG_ERROR, "Finger prints not retrieved\n");
        return;
    }
}

static int query_formats(AVFilterContext *ctx)
{
    static const enum AVSampleFormat sample_fmts[] = {
        AV_SAMPLE_FMT_DBL,
        AV_SAMPLE_FMT_DBLP,
        AV_SAMPLE_FMT_NONE
    };
    AVFilterFormats *formats;

    if (!(formats = ff_make_format_list(sample_fmts)))
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, formats);
}

static int filter_frame(AVFilterLink *inlink, AVFrame *samples)
{
    AVFilterContext *ctx = inlink->dst;
    PeakPointsContext *p = ctx->priv;
    int i, nb_samples = samples->nb_samples;

    // store audio data
    for (i = 0; i < nb_samples; i++) {
        p->data[p->index] = ((double*)(samples->data[0]))[i];
        p->buffFlag = 1;
        p->index = p->index + 1;

        // size check
        if (p->index == SIZECHECK) {
            // get peak points stats
            ppointsStats(ctx, p);
            p->time = p->time + 1;
            p->index = 0;
            p->buffFlag = 1;
        }
    }

    return ff_filter_frame(inlink->dst->outputs[0], samples);
}

static av_cold void uninit(AVFilterContext *ctx)
{
    PeakPointsContext *p = ctx->priv;

    // if audio data in buffer get peak points
    if (p->buffFlag) {
        ppointsStats(ctx, p);
    }

    // free allocated memories
    av_freep(&p->data);
    //av_freep(&p->peaks);
    av_freep(&p->fprints);
}

static const AVFilterPad peakpoints_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_AUDIO,
        .filter_frame = filter_frame,
    },
    { NULL }
};

static const AVFilterPad peakpoints_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_AUDIO,
    },
    { NULL }
};

AVFilter ff_af_peakpoints = {
    .name          = "peakpoints",
    .description   = NULL_IF_CONFIG_SMALL("peak points from frequency domain windowed data."),
    .init          = init,
    .uninit        = uninit,
    .query_formats = query_formats,
    .priv_size     = sizeof(PeakPointsContext),
    .inputs        = peakpoints_inputs,
    .outputs       = peakpoints_outputs,
    .priv_class    = &peakpoints_class,
};

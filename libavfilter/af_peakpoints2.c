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

#define SIZECHECK 4096*5
/* Upper and Lower Limit of sound frequency in audible range */
#define UPPER_LIMIT 300
#define LOWER_LIMIT 20
/* Total number of points in a single fingerprint, each point comes from a bin*/
#define KEYPOINTCOUNT 5


/* structure to store constellation point and time of its appearence */
typedef struct {
    size_t time;
    int frequency;
} ConstellationPoint;


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
    ConstellationPoint *cpoints;
    //double *peaks;
    int size; // number of peaks
    int windowSize;
    //char *inputFile;
} PeakPointsContext;


/* returns maximum value from an array conditioned on start and end index */
/*static double getMax(double *res_arr, int startIndex, int endIndex) {
    int i;
    double max = res_arr[startIndex];
    for (i = startIndex; i <= endIndex; i++) {
	    if (res_arr[i] > max) {
	        max = res_arr[i];
	    }
    }
    return max;
}*/

/* Stores peak frequency for each window(of chunkSize) in peaks array */
/*static void getPeakPointInChunk(int chunkSize, double *res_arr, int size, double *peaks) {
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
}*/

/* Find out the range in which the frequency point lies*/
static int getRangeIndex(AVFilterContext *ctx, int freq, int *range) {
    int i = 0;

    // null check
    /*if (!range) {
        av_log(ctx, AV_LOG_ERROR, "range is NULL\n");
    }*/

    while (range[i] < freq) {
        i++;
    }

    return i;
}

/* absolute value of complex number */
static double getAbs(FFTComplex cn) {
    //double res;

    //res = pow(cn.re, 2) + pow(cn.im, 2);
    //res = sqrt(res);
    return hypot(cn.re, cn.im);
}

/* Getting constellation points */
static ConstellationPoint *getConstellationPoints(AVFilterContext *ctx, PeakPointsContext *p,
        FFTComplex *tab, int *bins, int bin_size, size_t time) {
    int freq, index, j, start, end, flag, count;
    double mag, *maxscores;
    int *interval, *frequencies;
    ConstellationPoint *cpt = NULL;

    //av_log(ctx, AV_LOG_INFO, "Inside getFingerPrint method\n");

    // allocate memmory for cpt if not allocated

    cpt = av_malloc_array(bin_size, sizeof(*cpt));

    // memory check
    if (!cpt) {
        av_log(ctx, AV_LOG_ERROR, "Can't allocate memory for cpt array\n");
    }

    // initialize time and frequency to negative quantity
    index = 0;
    while (index < bin_size) {
        cpt[index].frequency = -1;
        cpt[index].time = -1;
        index = index + 1;
    }

    index = 0;

    interval = av_malloc_array(bin_size, sizeof(*interval));

    // memmory check
    if (!interval) {
        av_log(ctx, AV_LOG_ERROR, "Can't allocate memory for interval array\n");
    }
    //
    //av_log(ctx, AV_LOG_INFO, "copying\n");
    memcpy(interval, bins, bin_size * sizeof(int));
    //av_log(ctx, AV_LOG_INFO, "copying finished\n");

    maxscores = av_malloc_array(bin_size, sizeof(*maxscores));

    //set all elements to zero
    //av_log(ctx, AV_LOG_INFO, "setting\n");
    memset(maxscores, 0, bin_size * sizeof(*maxscores));

    /*for (int l =0; l < 5; l++) {
        av_log(ctx, AV_LOG_INFO, "%d, ", highscores[l]);
    }
    av_log(ctx, AV_LOG_INFO, "\n");*/

    //av_log(ctx, AV_LOG_INFO, "setting finished\n");

    // memmory check
    if (!maxscores) {
        av_log(ctx, AV_LOG_ERROR, "Can't allocate memmory for maxscore array\n");
    }

    frequencies = av_malloc_array(bin_size, sizeof(*frequencies));

    memset(frequencies, -1, bin_size * sizeof(*frequencies));

    // memory check
    if (!frequencies) {
        av_log(ctx, AV_LOG_ERROR, "Can't allocate memory for frequencies array\n");
    }

    //av_log(ctx, AV_LOG_INFO, "Before the loop\n");
    for (freq = 0; freq < p->windowSize; freq++) {
        // calculate magnitude in decibel
        //av_log(ctx, AV_LOG_INFO, "freq is %d\n", freq);
        mag = getAbs(tab[freq]);
        //av_log(ctx, AV_LOG_INFO, "got magnitude\n");

        // Find out the interval in which it lies
        index = getRangeIndex(ctx, freq, interval);
        //av_log(ctx, AV_LOG_INFO, "got index %d\n", index);

        // Save the highest magnitude and corresponding frequency
        if (mag > maxscores[index]) {
            maxscores[index] = mag;
            frequencies[index] = freq;
        }
        //av_log(ctx, AV_LOG_INFO, "saved\n");
    }

    // check for valid constellation points
    for (index = 0; index < bin_size; index++) {
        flag = 0;

        if (frequencies[index] < 64) {
            if (frequencies[index] == -1) {
                continue;
            }
            start = 0;
            end = frequencies[index] + 64;
        }
        else if (frequencies[index] > p->windowSize-64) {
            end = p->windowSize -1;
            start = frequencies[index] - 64;
        }
        else {
            start = frequencies[index] - 64;
            end = frequencies[index] + 64;
        }

        count = 1;
        for (j = start; j <= end; j++) {
            if (maxscores[index] < getAbs(tab[j])) {
                cpt[index].frequency = -1;
                cpt[index].time = -1;
                flag = 1;
            }

            if (maxscores[index] < getAbs(tab[j]) * 2) {
                count ++;
            }

            if (count > 3) {
                cpt[index].frequency = -1;
                cpt[index].time = -1;
                flag = 1;
            }
        }

        // case of valid points
        if (!flag) {
            cpt[index].frequency = frequencies[index];
            cpt[index].time = time;
        }
    }

    // free allocated memories.
    av_freep(&interval);
    av_freep(&maxscores);
    av_freep(&frequencies);
    return cpt;
}

/* Get peaks points from windowed frequency domain data*/
static int getPeakPoints2(AVFilterContext *ctx, PeakPointsContext *ppc) {
    int i, m, k, size, chunkSize, pSize, chunkSampleSize, resSize;
    int j, q, inverse, delta, bin_size, val, curr_index, t;
    int *bins;
    //double *fft_res;
    //void *avc;
    ConstellationPoint *cp;
    //RDFTContext *rdftC;
    FFTContext *fftc;
    FFTComplex *tab;
    //FFTSample *data;

    size = ppc->index;
    m = log2(ppc->windowSize);
    //fft_size = 1 << m;
    chunkSize = ppc->windowSize;
    chunkSampleSize = size/chunkSize;
    resSize = chunkSize * chunkSampleSize;
    pSize = resSize/chunkSize;

    // freqency bins
    delta = 64;
    val = 0;
    t = 0;

    // estimate bin size
    if ((chunkSize % delta) == 0) {
        bin_size = chunkSize/delta;
    } else {
        bin_size = chunkSize/delta + 1;
    }

    // allocating memory for bins array
    bins = av_malloc_array(bin_size, sizeof(*bins));

    // memory check
    if (!bins) {
        av_log(ctx, AV_LOG_ERROR, "Can't allocate memory for bins\n");
    }

    // fill bins with range values
    while (val < chunkSize) {
        bins[t] = val + delta;
        val =  val + delta;
        t = t + 1;
    }

    ppc->bin_size = bin_size;

    //fft_res = av_malloc_array(resSize, sizeof(double));

    /*if (!fft_res) {
        av_log(ctx, AV_LOG_ERROR, "Cann't allocate memmory for storing fft data\n");
        return 0;
    }*/

    ppc->cpoints = av_malloc_array(pSize*resSize, sizeof(ConstellationPoint));

    // memory check
    if (!ppc->cpoints) {
        av_log(ctx, AV_LOG_ERROR, "Can't allocate memory for cpoints\n");
    }

    //print chunk size
    //av_log(ctx, AV_LOG_INFO, "Chunk Size/Window Size is: %d", chunkSize);

    // fft_size or chunkSize
    tab = av_malloc(chunkSize * sizeof(FFTComplex));

    // memmory check
    if (!tab) {
        av_log(ctx, AV_LOG_ERROR, "Cann't allocate memory for complex data\n");
    }

    //rdftC = av_rdft_init(m, DFT_R2C);
    inverse = 0;
    fftc = av_fft_init(m, inverse);

    /*data = av_malloc_array(chunkSize, sizeof(FFTSample));

    if (!data) {
        av_log(ctx, AV_LOG_ERROR, "Cann't allocate memmory for chunk fft data\n");
        return 0;
    }*/

    // FFT transform for windowed time domain data
    // window is of size chunkSize
    k = 0;
    j = 0;
    curr_index = 0;
    //av_log(ctx, AV_LOG_INFO, "resSize is : %d\n", resSize);
    while (k < resSize) {
        //av_log(ctx, AV_LOG_INFO, "inside while loop\n");
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
        //av_log(ctx, AV_LOG_INFO, "calling getFingerprint method\n");
        cp = getConstellationPoints(ctx, ppc, tab, bins, bin_size, ppc->time);
        //av_log(ctx, AV_LOG_INFO, "after fingerprint method\n");
        ppc->time = ppc->time + 1;

        // assign constellation points to cpoints array
        q = 0;
        while (q < bin_size) {
            ppc->cpoints[curr_index + q] = cp[q];
            q = q + 1;
        }

        curr_index = curr_index + bin_size;
        //av_log(ctx, AV_LOG_INFO, "stored fp\n");

        j = j + 1;
        k = k + chunkSize;
        //av_log(ctx, AV_LOG_INFO, "j is %d and k is %d\n", j, k);
    }

    //av_log(ctx, AV_LOG_INFO, "loop end\n");

    //av_rdft_end(rdftC);
    av_fft_end(fftc);
    //pSize = resSize/chunkSize;
    // new size
    ppc->size = pSize * bin_size;
    //ppc->peaks = av_malloc_array(pSize, sizeof(double));

    /*if (!ppc->peaks) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memory for peak storage\n");
        return 0;
    }*/

    //getPeakPointInChunk(chunkSize, fft_res, resSize, ppc->peaks);

    // free allocated memories
    //av_freep(&data);
    av_freep(&tab);
    av_freep(&bins);
    //av_freep(&fft_res);
    //av_log(ctx, AV_LOG_INFO, "returning \n");
    return 1;
}

/* static void free_constellation_struct(AVFilterContext *ctx, ConstellationPoint *cpt, int size) {
    int i;

    // free keyPoints array first
    //av_log(ctx, AV_LOG_INFO, "freeing keypoints..\n");
    for (i = 0; i < size; i++) {
        av_freep(&cpt);
        free(&cpt[i]);
    }
    //av_log(ctx, AV_LOG_INFO, "freed keypoints\n");
    //av_freep(&fp);
    //av_log(ctx, AV_LOG_INFO, "freed fp\n");
}*/


#define OFFSET(x) offsetof(PeakPointsContext, x)

static const AVOption peakpoints2_options[] = {
    { "wsize",  "set window size", OFFSET(windowSize),  AV_OPT_TYPE_INT,    {.i64=1024},    0, INT_MAX},
    { NULL },
};

AVFILTER_DEFINE_CLASS(peakpoints2);

static av_cold int init(AVFilterContext *ctx)
{
    PeakPointsContext *p = ctx->priv;

    if (p->windowSize < 1024 || p->windowSize > SIZECHECK) {
        av_log(ctx, AV_LOG_ERROR, "window size must be in range 16 to %d\n", SIZECHECK);
        return AVERROR(EINVAL);
    }

    av_log(ctx, AV_LOG_INFO, "window size is %d\n", p->windowSize);

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
    ret = getPeakPoints2(ctx, p);

    if (ret && p->size) {
    	// print peaks
        if (p->isOnce) {
            av_log(ctx, AV_LOG_INFO, "######## Finger prints are ########\n");
            p->isOnce = 0;
        }

        for (i = 0; i < p->size; i++) {
            /*for (j = 0; j < p->bin_size; j++) {
                av_log(ctx, AV_LOG_INFO, "%d ", p->fprints[i].keyPoints[j]);
            }
            av_log(ctx, AV_LOG_INFO, "at timestep: %zu\n", p->fprints[i].time);
            */

            // print hashes for pair points
            // check for invalid points
            if (p->cpoints[i].frequency == -1) {
                continue;
            }

            for (j = i+1; j < p->size; j++) {
                // f1:f2:t2-t1:start time
                // check for invalid points
                if (p->cpoints[j].frequency == -1) {
                    continue;
                }

                av_log(ctx, AV_LOG_INFO, "%d:%d:%zu:%zu\n", p->cpoints[i].frequency,
                p->cpoints[j].frequency, abs(p->cpoints[j].time - p->cpoints[i].time),
                p->cpoints[i].time);
            }
        }

        // free peaks and set size to zero
        //av_freep(&p->peaks);
        av_freep(&p->cpoints); // check if memory to be freed hierarchial fashion
        //free_fingerprintstruct(ctx, p->cpoints, p->size);
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
            //p->time = p->time + 1;
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
    av_freep(&p->cpoints);
}

static const AVFilterPad peakpoints2_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_AUDIO,
        .filter_frame = filter_frame,
    },
    { NULL }
};

static const AVFilterPad peakpoints2_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_AUDIO,
    },
    { NULL }
};

AVFilter ff_af_peakpoints2 = {
    .name          = "peakpoints2",
    .description   = NULL_IF_CONFIG_SMALL("peak points from frequency domain windowed data."),
    .init          = init,
    .uninit        = uninit,
    .query_formats = query_formats,
    .priv_size     = sizeof(PeakPointsContext),
    .inputs        = peakpoints2_inputs,
    .outputs       = peakpoints2_outputs,
    .priv_class    = &peakpoints2_class,
};

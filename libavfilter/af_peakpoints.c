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

/* Structure to contain peak points context */
typedef struct {
    const AVClass *class;
    double *data;
    int nsamples;
    int index;
    int isOnce;
    int buffFlag;
    double *peaks;
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

/* Get peaks points from windowed frequency domain data*/
static int getPeakPoints(PeakPointsContext *ppc) {
    int i, m, k, size, chunkSize, pSize, chunkSampleSize, resSize;
    double *fft_res;
    void *avc;
    RDFTContext *rdftC;
    FFTSample *data;

    size = ppc->index;
    m = log2(ppc->windowSize);
    chunkSize = ppc->windowSize;
    chunkSampleSize = size/chunkSize;
    resSize = chunkSize * chunkSampleSize;

    fft_res = av_malloc_array(resSize, sizeof(double));

    if (!fft_res) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memmory for storing fft data\n");
        return 0;
    }


    rdftC = av_rdft_init(m, DFT_R2C);
    data = av_malloc_array(chunkSize, sizeof(FFTSample));

    if (!data) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memmory for chunk fft data\n");
        return 0;
    }
    // FFT transform for windowed time domain data
    // window is of size chunkSize
    k = 0;
    while (k < resSize) {
        //copy data
        for (i = 0; i < chunkSize; i++) {
            data[i] = ppc->data[i+k];
        }
        //calculate FFT
        av_rdft_calc(rdftC, data);
        for (i = 0; i < chunkSize; i++) {
            fft_res[i+k] = data[i];
        }
        k = k + chunkSize;
    }

    av_rdft_end(rdftC);
    pSize = resSize/chunkSize;
    ppc->size = pSize;
    ppc->peaks = av_malloc_array(pSize, sizeof(double));

    if (!ppc->peaks) {
        av_log(avc, AV_LOG_ERROR, "Cann't allocate memory for peak storage\n");
        return 0;
    }

    getPeakPointInChunk(chunkSize, fft_res, resSize, ppc->peaks);
    av_freep(&data);
    av_freep(&fft_res);
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
    int i, ret;
    ret = getPeakPoints(p);

    if (ret && p->size) {
    	// print peaks
        if (p->isOnce) {
            av_log(ctx, AV_LOG_INFO, "######## Peak points are ########\n");
            p->isOnce = 0;
        }
        for (i = 0; i < p->size; i++) {
            av_log(ctx, AV_LOG_INFO, "%f\n", p->peaks[i]);
        }
        // free peaks and set size to zero
        av_freep(&p->peaks);
        p->size = 0;
    } else if (p->size || !ret) {
        av_log(ctx, AV_LOG_ERROR, "Peak points not retrieved\n");
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
    av_freep(&p->peaks);
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

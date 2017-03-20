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
#include "libavcodec/bytestream.h"
#include "libavformat/avformat.h"
#include "libswscale/swscale.h"
#include "avfilter.h"
#include "audio.h"
#include "libavutil/opt.h"
#include "libavutil/internal.h"
#include "window_func.h"

#include <fcntl.h>

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

//#define SIZECHECK 4096*5
#define SIZECHECK 1024
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

/* structure to store match information */
typedef struct {
    int matchtime;
    int songid;
} MatchInfo;

/* Structure to contain peak points context */
typedef struct {
    const AVClass *class;
    double *data;
    int mode;
    int nsamples;
    int index;
    int song_id;
    size_t time;
    int isOnce;
    int buffFlag;
    int bin_size;
    float *window_func_lut;
    ConstellationPoint *cpoints;
    MatchInfo *mi;
    int mi_size;
    int mi_index;
    //double *peaks;
    int size; // number of peaks
    int windowSize;
    //char *inputFile;
} PeakPointsContext;


int cmpfunc(const MatchInfo *, const MatchInfo *);

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
        index++;
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
        else if (frequencies[index] >= p->windowSize-64) {
            end = p->windowSize -1;
            start = frequencies[index] - 64;
        }
        else {
            start = frequencies[index] - 64;
            end = frequencies[index] + 64;
        }

        count = 0;

        for (j = start; j <= end; j++) {
            if (maxscores[index] < getAbs(tab[j])) {
                cpt[index].frequency = -1;
                cpt[index].time = -1;
                flag = 1;
            }

            if (maxscores[index] < getAbs(tab[j]) * 2) {
                count ++;
            }

            if (count > 30) {
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
    int j, q, inverse, delta, bin_size, val, curr_index, t, lim;
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
        t++;
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
    i = 0;
    lim = 0;
    //av_log(ctx, AV_LOG_INFO, "resSize is : %d\n", resSize);
    while (k < resSize) {
        //av_log(ctx, AV_LOG_INFO, "inside while loop\n");
        //copy data; put imaginary part as 0
        for (i = 0; i < chunkSize; i++) {
            //data[i] = ppc->data[i+k];
            tab[i].re = ppc->data[i+lim] * ppc->window_func_lut[i];
            tab[i].im = 0.0;
        }

        // apply logic to get overlap ffts
        lim += 512;
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
        ppc->time++;

        // assign constellation points to cpoints array
        q = 0;
        while (q < bin_size) {
            ppc->cpoints[curr_index + q] = cp[q];
            q++;
        }

        curr_index = curr_index + bin_size;
        //av_log(ctx, AV_LOG_INFO, "stored fp\n");

        j++;
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
    { "mode", "do lookup for the song or store the song", OFFSET(mode), AV_OPT_TYPE_INT, {.i64=1}, 0, 1},
    { "song_id", "song identifier while storing a song", OFFSET(song_id), AV_OPT_TYPE_INT, {.i64=-1}, -1, INT_MAX},
    { NULL },
};

AVFILTER_DEFINE_CLASS(peakpoints2);

static av_cold int init(AVFilterContext *ctx)
{
    PeakPointsContext *p = ctx->priv;
    float overlap;

    if (p->windowSize < 1024 || p->windowSize > SIZECHECK) {
        av_log(ctx, AV_LOG_ERROR, "window size must be in range 16 to %d\n", SIZECHECK);
        return AVERROR(EINVAL);
    }

    if ((p->mode != 0) && (p->mode != 1)) {
        av_log(ctx, AV_LOG_ERROR, "mode is either 0 or 1; 0 to store, 1 to lookup\n");
    }

    // case of creating indexes
    if (p->mode == 0 && p->song_id == -1) {
        av_log(ctx, AV_LOG_ERROR, "Index creation requires a associated song_id\n");
    }

    av_log(ctx, AV_LOG_INFO, "window size is %d\n", p->windowSize);

    p->index = 0;
    p->time = 0;
    p->size = 0;
    p->isOnce = 1;
    p->mi_index = 0;
    p->mi = NULL;
    p->data = av_malloc_array(SIZECHECK, sizeof(double));

    p->window_func_lut = av_malloc_array(p->windowSize, sizeof(float));
    ff_generate_window_func(p->window_func_lut, p->windowSize , WFUNC_HAMMING, &overlap);

    if (!p->data) {
        av_log(ctx, AV_LOG_ERROR, "Cann't allocate memmory for audio data\n");
        return AVERROR(EINVAL);
    }

    return 0;
}

int cmpfunc(const MatchInfo * a, const MatchInfo * b) {
    int diff;

    diff = a->matchtime - b->matchtime;

    if (!diff) {
        return (a->songid - b->songid);
    }
    else {
        return diff;
    }
}

static void ppointsStats(AVFilterContext *ctx, PeakPointsContext *p) {
    int i, j, fp, ret, f1, f2, f3, f4, length, buf_size, song_id,
        val, count, retval, match_count, found, npeaks, songid, tstart;
    size_t t1, t2, t3, t4;
    char *filename;
    uint8_t entry[34];
    uint8_t buff[34];
    PutByteContext pb;
    GetByteContext gb;

    ret = getPeakPoints2(ctx, p);

    if (ret && p->size) {
    	// print peaks
        if (p->isOnce) {
            av_log(ctx, AV_LOG_INFO, "######## Finger prints are ########\n");
            p->isOnce = 0;
        }

        // get rid of invalid points
        /*mark_index = 0;
        for (i = 0; i < p->size; i++) {
            if (p->cpoints[i].frequency != -1) {
                    p->cpoints[mark_index] = p->cpoints[i];
                    mark_index = mark_index + 1;
            }

            else if (p->cpoints[i].frequency == -1) {
                mark_index = i;
            }
        }

        p->size = mark_index;*/

        /* allocate memory for MatchInfo array*/
        //p->mi = av_fast_realloc(p->mi, &p->mi_size, 1000);

        // check;
        /*if (!p->mi_size || !p->mi) {
            av_log(ctx, AV_LOG_ERROR, "MatchInfo array not reallocated\n");
        }*/

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

            /*for (j = i+1; j < p->size; j++) {
                // f1:f2:t2-t1:start time
                // check for invalid points
                if (p->cpoints[j].frequency == -1) {
                    continue;
                }

                av_log(ctx, AV_LOG_INFO, "%d:%d:%zu:%zu\n", p->cpoints[i].frequency,
                p->cpoints[j].frequency, abs(p->cpoints[j].time - p->cpoints[i].time),
                p->cpoints[i].time);
            }*/

            // considering 8 points now
            if (i+7 >= p->size) {
                break;
            }

            f1 = p->cpoints[i].frequency;
            t1 = p->cpoints[i].time;
            f2 = p->cpoints[i+1].frequency;
            t2 = p->cpoints[i+1].time;
            f3 = p->cpoints[i+2].frequency;
            t3 = p->cpoints[i+2].time;
            f4 = p->cpoints[i+3].frequency;
            t4 = p->cpoints[i+3].time;

            av_log(ctx, AV_LOG_INFO, "%d:%d:%d:%d:%zu:%zu:%zu\n", f1,
                    f2, f3, f4, t2-t1, t3-t2, t4-t3);

	        //filename is same as anchor point
            length = snprintf( NULL, 0, "%d", f1);
            length = length + snprintf(NULL, 0, "PPDB-%d", f1);
            filename = malloc(length + 1); //extra 1 for terminating character
            snprintf(filename, length + 1, "PPDB-%d", f1);

            if (!p->mode) {
                // put to file db
                buf_size = sizeof(entry);

                bytestream2_init_writer(&pb, entry, buf_size);
                // write data to putbyte context
                // skip entires that have less than 3 freq peaks
                npeaks = 0;
                for (j = 0; j < 8; j++) {
                    if (p->cpoints[i+j].frequency != -1) {
                        npeaks++;
                    }
                }

                // discard entries with less than 3 peak points
                if (npeaks < 2) {
                    continue;
                }

                for (j = 0; j < 8; j++) {
                    bytestream2_put_le16(&pb, p->cpoints[i+j].frequency);
                    /*if (j != 7) {
                        bytestream2_put_le16(pb, ',');
                    }*/
                }

                //bytestream2_put_le16(pb, ':');

                for (j = 0; j < 8; j++) {
                    bytestream2_put_le16(&pb, p->cpoints[i+j].time);
                    /*if (j != 7) {
                        bytestream2_put_le16(pb, ',');
                    }*/
                }

                // just for test song_id = 1
                song_id = p->song_id;
                //bytestream2_put_le16(pb, ':');
                bytestream2_put_le16(&pb, song_id);
                //bytestream2_put_le16(pb, '\n');

                fp = avpriv_open(filename, O_RDWR|O_CREAT|O_APPEND, 0666);

                //check
                if (fp == -1) {
                    av_log(ctx, AV_LOG_ERROR, "Cann't open file %s for writing\n", filename);
                }

                write(fp, entry, buf_size);

                close(fp);
                av_freep(&filename);
            }

            // lookup
            else if (p->mode) {
                // open corrosponding file and read the content in array
                fp = avpriv_open(filename, O_RDONLY);

                // file not found
                if (fp == -1) {
                    av_log(ctx, AV_LOG_ERROR, "No index file %s. Try storing it as index setting mode as 0\n", filename);
                    av_freep(&filename);
                    continue;
                }

                while (1) {
                    count = 0;
                    found = 0;
                    match_count = 0;

                    retval = read(fp, buff, sizeof(buff)); // is size to be read buf_size? should be total size of file content.

                    if (retval < 0) {
                        av_log(ctx, AV_LOG_ERROR, "No content read from file %s\n", filename);
                        //av_freep(&filename);
                        break;
                    }

                    if (retval == 0) {
                        //av_log(ctx, AV_LOG_INFO, "Reached EOF\n");
                        break;
                    }

                    bytestream2_init(&gb, buff, sizeof(buff));

                    // 17 is total number of items stored in one hash
                    for (count = 0; count < 8; count++) {
                        val = (int16_t)bytestream2_get_le16(&gb);
                        // match
                        /*if (!(match_time) && p->cpoints[i+count].frequency != val) {
                            // this hash doesn't match try next hash
                            break;
                        }
                        else if (!(match_time) && p->cpoints[i+count].frequency == val){
                            match_count = match_count + 1;
                        }

                        if (match_time) {
                            if (p->cpoints[i+count].time != val) {
                                break;
                            }
                            else {
                                match_count = match_count + 1;
                            }

                        }
                        if (match_count == 8) {
                            // all frequencies matched now try matching times
                            match_time = 1;
                            count = 0;
                        }
                        else if (match_count == 16) {
                            // it's a match
                            found = 1;
                            av_log(ctx, AV_LOG_INFO, "match found. ");
                        }*/
                        if ((val == -1) || (p->cpoints[i+count].frequency == -1)) {
                            continue;
                        }

                        if (p->cpoints[i+count].frequency != val) {
                            // this hash no a match, try next hash
                            //break;
                        }
                        else {
                            match_count = match_count + 1;
                        }

                        if ((match_count > 0)) {
                            //match found
                            //av_log(ctx, AV_LOG_INFO, "match found\n");
                            found = 1;
                        }
                    }

                    if (found) {
                        // skipping time values
                        for (count = 0; count < 8; count++) {
                            val = (int16_t)bytestream2_get_le16(&gb);

                            if (count == 0) {
                                tstart = val;
                            }
                        }

                        songid = (int16_t)bytestream2_get_le16(&gb);
                        //av_log(ctx, AV_LOG_INFO, "Song id is %d\n", songid);

                        p->mi = av_fast_realloc(p->mi, &p->mi_size, sizeof(MatchInfo)*(p->mi_index+1));

                        //check
                        if (!p->mi_size || !p->mi) {
                            av_log(ctx, AV_LOG_ERROR, "MatchInfo array not reallocated\n");
                        }
                        // store info in MatchInfo array
                        p->mi[p->mi_index].matchtime = p->cpoints[i].time - tstart;
                        p->mi[p->mi_index].songid = songid;
                        p->mi_index++;
                    }

                }

                close(fp);
                av_freep(&filename);
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
    int ret, sample_rates[] = { 11025, -1 };
    AVFilterChannelLayouts *layouts = NULL;

    static const enum AVSampleFormat sample_fmts[] = {
        AV_SAMPLE_FMT_DBL,
        AV_SAMPLE_FMT_DBLP,
        AV_SAMPLE_FMT_NONE
    };
    AVFilterFormats *formats;

    formats = ff_make_format_list(sample_rates);

    if (!formats)
        return AVERROR(ENOMEM);

    ret = ff_set_common_samplerates(ctx, formats);

    if (ret < 0)
        return ret;

    ret = ff_add_channel_layout(&layouts, AV_CH_LAYOUT_MONO);

    if (ret)
        return ret;

    if (!layouts)
        return AVERROR(ENOMEM);

    ret = ff_set_common_channel_layouts(ctx, layouts);

    if (ret < 0)
        return ret;

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
        p->index++;

        // size check
        if (p->index == SIZECHECK) {
            // get peak points stats
            ppointsStats(ctx, p);
            //p->time = p->time + 1;
            memmove(p->data, p->data + p->windowSize/2, sizeof(p->data)*SIZECHECK/2);
            p->index = 512;
            p->buffFlag = 1;
        }
    }

    return ff_filter_frame(inlink->dst->outputs[0], samples);
}

static av_cold void uninit(AVFilterContext *ctx)
{
    int i, curr_max, count, prev_songid, curr_songid, match_songid;
    int prev_matchtime, curr_matchtime;
    PeakPointsContext *p = ctx->priv;

    // if audio data in buffer get peak points
    if (p->buffFlag) {
        ppointsStats(ctx, p);
    }

    // in lookup mode and any match found
    if (p->mode && p->mi_index) {
        // sort matchinfo array based on timediff
        qsort(p->mi, p->mi_index, sizeof(MatchInfo), cmpfunc);

        // find longest consecutive run of identical songid
        count = 0;
        curr_max = 0;

        /*if (p->mi) {
            av_log(ctx, AV_LOG_INFO, "mi has content\n");
        }
        else { av_log(ctx, AV_LOG_INFO, "mi is null\n");}*/

        prev_songid = p->mi[0].songid;
        prev_matchtime = p->mi[0].matchtime;

        for (i = 0; i < p->mi_index; i++) {
            curr_songid = p->mi[i].songid;
            curr_matchtime = p->mi[i].matchtime;

            if ((prev_songid == curr_songid) && (prev_matchtime == curr_matchtime)) {
                count++;
                if (i == p->mi_index-1) {
                    if (count > curr_max) {
                        match_songid = prev_songid;
                        curr_max = count;
                    }
                }
            }
            else {
                if (count > curr_max) {
                    match_songid = prev_songid;
                    curr_max = count;
                }
                count = 1;
            }

            prev_songid = curr_songid;
            prev_matchtime = curr_matchtime;
        }

        // print match songid
        av_log(ctx, AV_LOG_INFO, "Matched Song ID is %d\n", match_songid);

        // free
        av_freep(&p->mi);
    }
    else if(p->mode && (!p->mi_index)) {
        av_log(ctx, AV_LOG_INFO, "No match found\n");
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

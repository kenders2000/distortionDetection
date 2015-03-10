/*Copyright (c) <2014> <Paul Kendrick>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */


#include "processInWindows.hpp"

struct wavfile {
    char id[4]; // should always contain "RIFF"
    int32_t totallength; // total file length minus 8
    char type[4]; // should always contain "WAVE "
    char typeId[4]; // should always contain "fmt "
    int32_t fmtchunksize; // format chunk size
    int16_t compression; // 1 pcm 3 ieee
    int16_t nochan; //
    int32_t fs; // fs
    int32_t bytes_per_sec; // byes / sec
    int16_t blockalign; //
    int16_t bits_per_sample; //
} __attribute__((__packed__));

static struct wavfile header;

void wavRms(char *filename, float *maxL, float* rms)
{
 FILE * wav;
  float tmp3;            


    wav = fopen(filename, "r");
    //check openable
    int test = sizeof (header);
    if (fread(&header, sizeof (header), 1, wav) < 1) {
        fprintf(stderr, "Can't read input file header %s\n", filename);
        exit(1);
    }
    printf("\nWavefile Header");
    printf("chunk ID %.*s\n", 4, header.id);
    printf("Total length %i\n", header.totallength);
    printf("should say WAVE %.*s\n", 4, header.type);
    printf("should say fmt : %.*s\n", 4, header.typeId);
    printf("format chunk size : %i\n", header.fmtchunksize);
    printf("Compression code %i\n", header.compression);
    printf("No. Channels %i\n", header.nochan);
    printf("Fs  %i \n", header.fs);
    printf("bytes_per_sec  %i \n", header.bytes_per_sec);
    printf("blockalign  %i\n", header.blockalign);
    printf("bits_per_sample  %i\n", header.bits_per_sample);

    if (header.fs != 44100.0) {
        fprintf(stderr, "Input must be sampled at 44100 Hz \n");
        exit(1);
    }
    // Skip though rest of header to find data 
    char datahead[4];
    int32_t datasize;
    int skip = 0;
    while (skip == 0) {
        fread(&datahead, sizeof (datahead), 1, wav);
        //  Skips through unitll it finds 'data'
        skip = memcmp(datahead, "data", 4) == 0;
        if (skip == 0) {
            fseek(wav, -3, SEEK_CUR);
            //printf("should say data %.*s\n", 4, datahead);
        }
    }

    // read data header 
    // printf("should say data %.*s\n", 4, datahead);
    fread(&datasize, sizeof (int32_t), 1, wav);
    printf(" datasize %i\n", datasize);            


    int32_t NoBlocks = datasize / (header.bits_per_sample / 8 * header.nochan);
            int i, c, n;   

            int frame=0;
        for (i = 0; i < (NoBlocks); i++) {
            //printf("%f ",(float)i/(float)NoBlocks);
            //channelsum[i]=0;
            for (c = 0; c < (header.nochan); c++) {
                float tmp2;
                float normV;
                if (header.bits_per_sample == 16 && header.compression == 1) {
                    int16_t tmp;
                    fread(&tmp, sizeof (int16_t), 1, wav);
                    tmp2 = (float) tmp;
                    normV = (65535.0 / 2.0);
                } else if (header.bits_per_sample == 32 && header.compression == 1) {
                    int32_t tmp = 0;
                    fread(&tmp, sizeof (int32_t), 1, wav);
                    tmp2 = (float) tmp;
                    normV = (2147483648.0);
                } else if (header.bits_per_sample == 32 && header.compression == 3) {
                    float tmp = 0;
                    fread(&tmp, sizeof (float), 1, wav);
                    tmp2 = tmp;
                    normV = 1; //(16777216.0/2.0);
                } else {
                    printf("Supported formats are 16 and 32 bit signed integer, and IEEE float little-endian - any number of channels");
                    break;
                }                     
                 tmp3 = tmp3 + pow(abs(tmp2 / (float) header.nochan),2); //*  158489=  10^(104/20)
                 if(abs(tmp3)> maxL[0])
                     maxL[0]=abs(tmp3);


                frame++;
            }
}

            rms[0]=sqrtf(tmp3/(double(frame))); 

}

void loadWav(char * filename, char * outFilename, const char *jsonFilename, char *treeDir, float gain, int frameAve, float thresh) {
    char str1[100], str2[100];
    FILE * pFile;
    FILE *pFileJSON;


//    sprintf(str1, "trees/%s/levelClass", treeDir);
//    sprintf(str2, "trees/%s/snrClass", treeDir);
    sprintf(str1, "trees/dectrees_new/trainedTree_96_bags");
    DTree distTree;
    printf("\nLoading Decision Trees");
    distTree.readTextFilesTrees(str1);
    printf("\nDone");

    //http://www.sonicspot.com/guide/wavefiles.html
    //http://yannesposito.com/Scratch/en/blog/2010-10-14-Fun-with-wav/char 
    //char  filename[]="test16S.wav";
    //char  filename[]="iphone1.wav";
    float* maxL=(float*) malloc(sizeof (float)*1);maxL[0]=1.5;                 

    float *rms=(float*) malloc(sizeof (float)*1);rms[0]=0;    
    wavRms(filename,maxL,rms);
    FILE * wav;
    wav = fopen(filename, "r");
    //check openable
    int test = sizeof (header);
    if (fread(&header, sizeof (header), 1, wav) < 1) {
        fprintf(stderr, "Can't read input file header %s\n", filename);
        exit(1);
    }
    printf("\nWavefile Header");
    printf("chunk ID %.*s\n", 4, header.id);
    printf("Total length %i\n", header.totallength);
    printf("should say WAVE %.*s\n", 4, header.type);
    printf("should say fmt : %.*s\n", 4, header.typeId);
    printf("format chunk size : %i\n", header.fmtchunksize);
    printf("Compression code %i\n", header.compression);
    printf("No. Channels %i\n", header.nochan);
    printf("Fs  %i \n", header.fs);
    printf("bytes_per_sec  %i \n", header.bytes_per_sec);
    printf("blockalign  %i\n", header.blockalign);
    printf("bits_per_sample  %i\n", header.bits_per_sample);

    if (header.fs != 44100.0) {
        fprintf(stderr, "Input must be sampled at 44100 Hz \n");
        exit(1);
    }
    // Skip though rest of header to find data 
    char datahead[4];
    int32_t datasize;
    int skip = 0;
    while (skip == 0) {
        fread(&datahead, sizeof (datahead), 1, wav);
        //  Skips through unitll it finds 'data'
        skip = memcmp(datahead, "data", 4) == 0;
        if (skip == 0) {
            fseek(wav, -3, SEEK_CUR);
            //printf("should say data %.*s\n", 4, datahead);
        }
    }

    // read data header 
    // printf("should say data %.*s\n", 4, datahead);
    fread(&datasize, sizeof (int32_t), 1, wav);

    printf("datasize %i\n", datasize);
    int32_t NoBlocks = datasize / (header.bits_per_sample / 8 * header.nochan);

    // allocate memory for time domain items
    double *windowD, *windowOverD, *windowPrevD;
    float *window;
    float *windowOver;
    float *windowPrev;
    float *last1s;
    float global_scale_max_over_1s = 0;
    float lastscale = 0;
    last1s = (float*) malloc(sizeof (float)*WIN_N * frameAve);
    int i;
    for (i = 0; i < (WIN_N * frameAve); i++)
        last1s[i] = 0;

    window = (float*) malloc(sizeof (float)*WIN_N);

    windowOver = (float*) malloc(sizeof (float)*WIN_N);
    windowPrev = (float*) malloc(sizeof (float)*WIN_N);

    windowD = (double*) malloc(sizeof (double)*WIN_N);
    windowPrevD = (double*) malloc(sizeof (double)*WIN_N);
    windowOverD = (double*) malloc(sizeof (double)*WIN_N);

    float * mfcc = (float*) malloc((16) * sizeof (float));
    float * mfccO = (float*) malloc((16) * sizeof (float));
    float * mfccP = (float*) malloc((16) * sizeof (float));

    double *feats;
    feats = (double*) malloc(sizeof (double)*800);
    double *featstmp;
    featstmp = (double*) malloc(sizeof (double)*800);


    int wN = 0;
    int N1s = 0;
    int wNo = WIN_N / 2;

    //temp features
    double *outMatrix = (double*) calloc(WIN_N, sizeof (double));
    double *lastspectrum = (double*) calloc(WIN_N, sizeof (double));
    double *spectrum = (double*) calloc(WIN_N, sizeof (double));
    double *lastOverspectrum = (double*) calloc(WIN_N, sizeof (double));
    int count_[1];
    count_[0] = 0;
    double *centroid = (double*) calloc(1, sizeof (double));
    double *spread = (double*) calloc(1, sizeof (double));
    double *kurt = (double*) calloc(1, sizeof (double));
    double *ent = (double*) calloc(1, sizeof (double));
    double *roughness = (double*) calloc(1, sizeof (double));
    double *skewness = (double*) calloc(1, sizeof (double));
    double*peaks = (double*) malloc(sizeof (double)*WIN_N);
    //mean features
    double *outMatrix_ = (double*) calloc(WIN_N, sizeof (double));
    int count__[1];
    count__[0] = 0;
    double *centroid_ = (double*) calloc(1, sizeof (double));
    double *spread_ = (double*) calloc(1, sizeof (double));
    double *kurt_ = (double*) calloc(1, sizeof (double));
    double *ent_ = (double*) calloc(1, sizeof (double));
    double *roughness_ = (double*) calloc(1, sizeof (double));
    double *skewness_ = (double*) calloc(1, sizeof (double));
    double*peaks_ = (double*) malloc(sizeof (double)*WIN_N);
    double* SpecFluxV = (double*) malloc(WIN_N * sizeof (double));
    double SpecFlux = 0;
    double SpecFlux2 = 0;
    sprintf(str2, "%s.b", outFilename);
    remove(str2);



    pFile = fopen(str2, "w");
    FILE * pFileRawData;
    //     char ss[100];
    //     sprintf(ss,"rawdata");     
    remove("rawdata.txt");
    pFileRawData = fopen("rawdata.txt", "w");
    FILE *pFileMeanData;
    remove("meandata.txt");
    pFileMeanData = fopen("meandata.txt", "w");
    //if the wav is 16 bit PCM
    int counter;
    int firstWin = 0;
    int classLevel;
    int classSNR;
    int frame = -1;
    int sumclassSNR = 0;
    int sumclassLevel = 0;
    float avermsAW = 0;
    initVarsFeats(WIN_N, header.fs);
    for (i = 0; i < (WIN_N); i++)
        window[i] = 0;

    if (1 == 1)//(header.bits_per_sample==16 && header.compression==1)
    {



        int i, c, n;
        for (i = 0; i < (NoBlocks); i++) {
            //printf("%f ",(float)i/(float)NoBlocks);
            //channelsum[i]=0;
            for (c = 0; c < (header.nochan); c++) {
                float tmp2;
                float normV;
                if (header.bits_per_sample == 16 && header.compression == 1) {
                    int16_t tmp;
                    fread(&tmp, sizeof (int16_t), 1, wav);
                    tmp2 = (float) tmp;
                    normV = (65535.0 / 2.0);
                } else if (header.bits_per_sample == 32 && header.compression == 1) {
                    int32_t tmp = 0;
                    fread(&tmp, sizeof (int32_t), 1, wav);
                    tmp2 = (float) tmp;
                    normV = (2147483648.0);
                } else if (header.bits_per_sample == 32 && header.compression == 3) {
                    float tmp = 0;
                    fread(&tmp, sizeof (float), 1, wav);
                    tmp2 = tmp;
                    normV = 1; //(16777216.0/2.0);
                } else {
                    printf("Supported formats are 16 and 32 bit signed integer, and IEEE float little-endian - any number of channels");
                    break;
                }


                float tmp3 = (tmp2 / (float) header.nochan); //*  158489=  10^(104/20)

                window[wN] = last1s[N1s]; //last window before overwritten
                last1s[N1s] = tmp3; // Overwrite

                if (abs(last1s[N1s]) > global_scale_max_over_1s)
                    global_scale_max_over_1s = abs(last1s[N1s]);
                // printf("%f \n",window[wN]);
                N1s++;
                wN++;

                if (N1s == (frameAve * WIN_N)) {
                    N1s = 0;
                    lastscale = global_scale_max_over_1s;
                    global_scale_max_over_1s = 0;
                }
            }

            //every WIN_N samples do this
            if (wN == WIN_N) {
                frame++;
                double scaleover_temp = 0;
                double scaleprev_temp = 0;
                double scale_temp = 0;
                double global_scale = 1;
                for (n = 0; n < (WIN_N / 2); n++) {
                    windowOver[n] = windowPrev[n + WIN_N / 2];
                }
                for (n = (WIN_N / 2); n < WIN_N; n++) {
                    windowOver[n] = window[n - WIN_N / 2];
                }
                //Normalise
                for (n = 0; n < WIN_N; n++) {
                    scaleover_temp = scaleover_temp + pow((double) windowOver[n], 2) / WIN_N;
                    scaleprev_temp = scaleprev_temp + pow((double) windowPrev[n], 2) / WIN_N;
                    scale_temp = scale_temp + pow((double) window[n], 2) / WIN_N;
                }
                scaleover_temp = sqrt(scaleover_temp);
                scaleprev_temp = sqrt(scaleprev_temp);
                scale_temp = sqrt(scale_temp);
                for (n = 0; n < WIN_N; n++) {
                    
//                tdData=tdData./sqrt(mean((tdData.^2)));
//                scale=max(abs(x1s(seci,:)./sqrt(mean(x1s(seci,:).^2))))/rms_store(i);
                    windowOverD[n] = ((double) windowOver[n] ) / scaleover_temp;
                    windowPrevD[n] = ((double) windowPrev[n] ) / scaleprev_temp;
                    windowD[n] = ((double) window[n] ) / scale_temp;
                }

                for (n = 0; n < WIN_N; n++) {
                    peaks[n] = 0;
                }
                // Previous window
                for (n = 0; n < 255; n++)
                    outMatrix[n] = 0;
                
                //scaling for the histogram is required using the max level over
                // the 1 s win which is lastscale, and we need to un-normalise the
                //  1024 sample window (was normed to rms level over 1024)
                gatewayFunctionDist(windowOverD,
                        WIN_N,
                        (double) header.fs,
                        lastscale / scaleover_temp,
                        outMatrix,
                        lastOverspectrum,
                        count_,
                        centroid,
                        spread,
                        kurt,
                        ent,
                        peaks,
                        roughness,
                        skewness);

                for (n = 0; n < 255; n++) {
                    featstmp[n] = outMatrix[n];
                }
                featstmp[255] = skewness[0];
                featstmp[256] = kurt[0];
                featstmp[257] = ent[0];
                featstmp[258] = roughness[0];
                featstmp[259] = count_[0];
                featstmp[260] = (SpecFlux2);

                for (n = 0; n < 261; n++) {
                    fprintf(pFileRawData, "%f ", featstmp[n]);
                }
                fprintf(pFileRawData, "\n");



                for (n = 0; n < 255; n++)
                    outMatrix_[n] = outMatrix_[n] + outMatrix[n];
                centroid_[0] = centroid_[0] + centroid[0] / ((double) frameAve * 2.0);
                spread_[0] = spread_[0] + spread[0] / ((double) frameAve * 2.0);
                kurt_[0] = kurt_[0] + kurt[0] / ((double) frameAve * 2.0);
                ent_[0] = ent_[0] + ent[0] / ((double) frameAve * 2.0);
                roughness_[0] = roughness_[0] + roughness[0] / ((double) frameAve * 2.0);
                skewness_[0] = skewness_[0] + skewness[0] / ((double) frameAve * 2.0);
                count__[0] = count__[0]+(count_[0] / (double) ((double) WIN_N / (double) header.fs * 2)) / (frameAve * 2);
                ;
                ;


                for (n = 0; n < WIN_N; n++) {
                    peaks[n] = 0;
                }
                //Overlapping window
                for (n = 0; n < 255; n++)
                    outMatrix[n] = 0;
                gatewayFunctionDist(windowD,
                        WIN_N,
                        (double) header.fs,
                        lastscale / scale_temp,
                        outMatrix,
                        spectrum,
                        count_,
                        centroid,
                        spread,
                        kurt,
                        ent,
                        peaks,
                        roughness,
                        skewness);

                for (n = 0; n < 255; n++) {
                    featstmp[n] = outMatrix[n];
                }
                featstmp[255] = skewness[0];
                featstmp[256] = kurt[0];
                featstmp[257] = ent[0];
                featstmp[258] = roughness[0];
                featstmp[259] = count_[0];
                featstmp[260] = (SpecFlux2);

                for (n = 0; n < 261; n++) {
                    fprintf(pFileRawData, "%f ", featstmp[n]);
                }
                fprintf(pFileRawData, "\n");

                for (n = 0; n < 255; n++)
                    outMatrix_[n] = outMatrix_[n] + outMatrix[n];
                centroid_[0] = centroid_[0] + centroid[0] / ((double) frameAve * 2.0);
                spread_[0] = spread_[0] + spread[0] / ((double) frameAve * 2.0);
                kurt_[0] = kurt_[0] + kurt[0] / ((double) frameAve * 2.0);
                ent_[0] = ent_[0] + ent[0] / ((double) frameAve * 2.0);
                roughness_[0] = roughness_[0] + roughness[0] / ((double) frameAve * 2.0);
                skewness_[0] = skewness_[0] + skewness[0] / ((double) frameAve * 2.0);
                count__[0] = count__[0]+(count_[0] / (double) ((double) WIN_N / (double) header.fs * 2)) / (frameAve * 2);



                SpecFlux = 0;

                for (n = 0; n < WIN_N / 2 + 1; n++)
                    SpecFlux = SpecFlux + ((pow(spectrum[n] - lastOverspectrum[n], 2)));
                SpecFlux2 = SpecFlux2 + sqrt(SpecFlux) / ((double)frameAve * 2.0);
                //            printf("%f \n",sqrt(SpecFlux));
                SpecFlux = 0;
                for (n = 0; n < WIN_N / 2 + 1; n++)
                    SpecFlux = SpecFlux + ((pow(lastOverspectrum[n] - lastspectrum[n], 2)));
                SpecFlux2 = SpecFlux2 + sqrt(SpecFlux) / ((double)frameAve * 2.0);
                counter++;
                if ((frame % (frameAve)) == 0 && (frame>=frameAve*2) ) {
                    double pmfmean=0;
                    for (n = 0; n < 255; n++) {
                        pmfmean=pmfmean+outMatrix_[n];
                    }

                    for (n = 0; n < 255; n++) {
                        feats[n] = outMatrix_[n]/pmfmean;
                        //printf("%f \n",feats[n]);
                    }
                    feats[255] = skewness_[0];
                    feats[256] = kurt_[0];
                    feats[257] = ent_[0];
                    feats[258] = roughness_[0];
                    feats[259] = count__[0];
                    feats[260] = (SpecFlux2);
                    //                printf("%f \n",feats[260]);
                    for (n = 0; n < 261; n++) {
                        fprintf(pFileMeanData, "%f ", feats[n]);
                        //printf("%f ",feats[n]);
                    }
                    int distLevel = distTree.decisionTreeFun(feats);
                    printf("frame centre %f Distortion Amount %i \n",(frame-frameAve/2.0),5-distLevel);
                    fprintf(pFileMeanData, "\n");
                    counter = 0;
                    //reset counters /vectors

                    for (n = 0; n < 255; n++)
                        outMatrix_[n] = 0;

                    centroid_[0] = 0;
                    spread_[0] = 0;
                    kurt_[0] = 0;
                    ent_[0] = 0;
                    roughness_[0] = 0;
                    skewness_[0] = 0;
                    count__[0] = 0;
                    SpecFlux2 = 0;

                }

                for (n = 0; n < WIN_N; n++) {
                    lastspectrum[n] = spectrum[n];
                }
                for (n = 0; n < WIN_N; n++)
                    windowPrev[n] = window[n];
                wN = 0;
      

            }



            // printf("%f \n",channelsum[i]);
            //fprintf(pFile,"%f \n",channelsum[i]);

        }
    }
    fclose(pFile);
    fclose(pFileRawData);
    fclose(pFileMeanData);

    // Global stats
    float aveAveLevel = 0;
    float maxAveLevel = 0;
    float minAveSNR = 5;
    float aveAveSNR = 0;
    //float thresh = 0.2;
    float avecomb = 0;
    float maxcomb = 0;
    float count_wF5 = 0;
    float count_wF4 = 0;
    float count_wF3 = 0;
    float count_wF2 = 0;
    float count_wF1 = 0;
    float count_wF0 = 0;
    counter = 0;
    //  pFile2 = fopen ( 'tmp.txt' , "w");
    pFile = fopen(str2, "r");
    char mystring[1000];
    while (fgets(mystring, sizeof (mystring), pFile) != NULL) {

        float dBA, t, level, snr, comb;
        sscanf(mystring, "%f %f %f %f %f", &t, &dBA,
                &level, &snr, &comb);

        if (level > maxAveLevel)
            maxAveLevel = level;

        if (snr < minAveSNR)
            minAveSNR = snr;

        if (comb > maxcomb)
            maxcomb = comb;

        avecomb = avecomb + comb;
        aveAveSNR = aveAveSNR + snr;
        aveAveLevel = aveAveLevel + level;

        if (level < 1)
            count_wF0++;


        if (snr >= 0 && snr < 1 && level >= 1)
            count_wF5++;
        if (snr >= 1 && snr < 2 && level >= 1)
            count_wF4++;
        if (snr >= 2 && snr < 3 && level >= 1)
            count_wF3++;
        if (snr >= 3 && snr < 4 && level >= 1)
            count_wF2++;
        if (snr >= 4 && snr < 5 && level >= 1)
            count_wF1++;



        //printf("%f %f %f %f %f \n", t, dBA, 
        //       level,snr, comb);
        counter++;
    }
    fclose(pFile);
    count_wF0 = count_wF0 / counter * 100.0;
    count_wF1 = count_wF1 / counter * 100.0;
    count_wF2 = count_wF2 / counter * 100.0;
    count_wF3 = count_wF3 / counter * 100.0;
    count_wF4 = count_wF4 / counter * 100.0;
    count_wF5 = count_wF5 / counter * 100.0;

    avecomb = avecomb / counter;
    aveAveSNR = aveAveSNR / counter;
    aveAveLevel = aveAveLevel / counter;

    // remove(outFilename);
    pFileJSON = fopen(jsonFilename, "w");
    //
    FILE * pFile2;
    pFile2 = fopen(outFilename, "w");
    pFile = fopen(str2, "r");
    fprintf(pFile2, "Microphone Wind Noise Detection - University of Salford - the Good Recording Project http://www.goodrecording.net \n\n", filename);
    fprintf(pFile2, "Microphone Wind noise Analysis for input file %s\n\n", filename);
    fprintf(pFile2, "Wind Noise Statistics, %% number of frames with wind noise detected at the following Signal to noise Ratios (high values = good quality low values = bad quality). also presented is the range in Quality Degradation that each SNR range represents.\n\n");
    fprintf(pFile2, "SNR range;\t>20dB,\t20 to 10dB,\t10 to 0dB,\t0 to -10dB,\t-10 to -20dB,\t<-20dB\n", count_wF0, count_wF1, count_wF2, count_wF3, count_wF4, count_wF5);
    fprintf(pFile2, "Qual. Degradation;\t<24%%,\t24 to 40%%,\t40 to 54%%,\t54 to 69%%,\t69 to 95%%,\t>95%%\n");
    fprintf(pFile2, "%% of time in each SNR range;\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f\n", count_wF0, count_wF1, count_wF2, count_wF3, count_wF4, count_wF5);
    fprintf(pFile2, "\n");
    fprintf(pFile2, "Wind noise time history\n\n");
    fprintf(pFile2, "T(s)\t\tQuality Degradation(%%)\t dBA \n");

    //JSON global stats
    fprintf(pFileJSON, "[\n\t{\n\t\"test\": \"Global Stats\",\n\t\"results\": [");
    fprintf(pFileJSON, "%0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f]\n\t},", count_wF0, count_wF1, count_wF2, count_wF3, count_wF4, count_wF5);

    fprintf(pFileJSON, "\n\t{\n\t\"test\": \"Time History\",\n\t\"results\": [\n");

    while (fgets(mystring, sizeof (mystring), pFile) != NULL) {

        float dBA, t, level, snr, comb;
        sscanf(mystring, "%f %f %f %f %f;", &t, &dBA,
                &level, &snr, &comb);
        float x = snr;
        float qual = 0.0093616465 * pow(x, 9) + -0.2408574076 * pow(x, 8) + 2.6409491876 * pow(x, 7) + -16.0598328947 * pow(x, 6) + 59.1080787296 * pow(x, 5) + -135.0673590478 * pow(x, 4) + 189.6816801547 * pow(x, 3) + -156.8624014598 * pow(x, 2) + 83.7429623557 * pow(x, 1) + 4.5226369050;
        if (comb == 0)
            qual = 100;

        fprintf(pFile2, "%4.1f\t\t%0.0f\t\t\t\t\t\t%4.0f\n", t, 100 - qual, dBA);
        if (firstWin == 1)
            fprintf(pFileJSON, ",\n");
        fprintf(pFileJSON, "\t\t{\"T\": %0.2f, \"dBA\": %0.0f, \"QDeg\": %0.2f}", t, dBA, 100 - qual);
        firstWin = 1;
        counter++;
    }

    ///////////////////////////// Find contiguous regions without wind noise
    fclose(pFile);
    pFile = fopen(str2, "r");
    fprintf(pFileJSON, "\n \t\t] \n \t},\n");
    fprintf(pFileJSON, "\t{\n\t\"test\": \"wind free regions\",\n\t\"results\": [\n");




    firstWin = 0;
    int clean = 1; // 
    int initialised = 0;
    int first = 0;
    float win = 0;
    float start = 0;
    fprintf(pFile2, "\nWind free regions from - to (s) using a Threshold of %2.0f\n\n", thresh);

    while (fgets(mystring, sizeof (mystring), pFile) != NULL) {

        float dBA, t, level, snr, comb;

        sscanf(mystring, "%f %f %f %f %f", &t, &dBA,
                &level, &snr, &comb);
        float x = snr;
        float qual = 0.0093616465 * pow(x, 9) + -0.2408574076 * pow(x, 8) + 2.6409491876 * pow(x, 7) + -16.0598328947 * pow(x, 6) + 59.1080787296 * pow(x, 5) + -135.0673590478 * pow(x, 4) + 189.6816801547 * pow(x, 3) + -156.8624014598 * pow(x, 2) + 83.7429623557 * pow(x, 1) + 4.5226369050;
        if (comb == 0)
            qual = 100;
        comb = 100 - qual;

        if (comb >= thresh && clean == 1 && initialised == 1 && first == 1) // previous frame was clean now its windy
        {
            fprintf(pFile2, "%0.2f\t%0.2f\n", start, t - win * 2);
            if (firstWin == 1)
                fprintf(pFileJSON, ",\n");
            fprintf(pFileJSON, "\t\t{\"s\": %0.2f, \"e\": %0.0f}", t, start, t - win * 2);
            firstWin = 1;

            clean = 0;

        } else if (comb < thresh && clean == 0 && initialised == 1) // previous frame was windy now its clean
        {
            start = t - win * 2;
            clean = 1;
            first = 1;
        }

        if (initialised == 0) {
            win = t / 2;

            if (comb >= thresh) {
                clean = 0;

            } else {
                clean = 1;
                first = 1;
                start = t - win * 2;
            }
            initialised = 1;
        }
        counter++;
    }
    ///////////////////////////

    fclose(pFile);
    fclose(pFile2);
    fprintf(pFileJSON, "\t]\n}\n]");
    fclose(pFileJSON);

    //remove(str2);

}


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
  //  int32_t NoBlocks;
} __attribute__((__packed__));
int32_t NoBlocks;
static struct wavfile header;

void wavRms(char *filename, float *maxL, float* rms) {
    FILE * wav;
    float tmp3;
    int verbose=0;

    wav = fopen(filename, "r");
     if(wav==NULL)
     {
    if (verbose ==1){   printf("Can't read input file \n");}
        exit(1);
    }
       
    //check openable
    int test = sizeof (header);
    if (fread(&header, sizeof (header), 1, wav) < 1) {
        fprintf(stderr, "Can't read input file header %s\n", filename);
       if (verbose ==1){  printf("Can't read input file header %s\n", filename);}

        exit(1);
    }
    if (verbose ==1){
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
    }
    if (header.fs != 44100.0) {
        fprintf(stderr, "Input must be sampled at 44100 Hz \n");
        exit(0);
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
    //printf("should say data %.*s\n", 4, datahead);
    // read data header 
    fread(&datasize, sizeof (int32_t), 1, wav);
//    printf(" datasize %i\n", datasize);


     NoBlocks = datasize / (header.bits_per_sample / 8 * header.nochan);
    //header.NoBlocks=NoBlocks;
    int i, c, n;

    int frame = 0;
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
//                printf("Supported formats are 16 and 32 bit signed integer, and IEEE float little-endian - any number of channels");
                break;
            }
            tmp3 = tmp3 + pow(abs(tmp2 / (float) header.nochan), 2); //*  158489=  10^(104/20)
            if (abs(tmp3) > maxL[0])
                maxL[0] = abs(tmp3);


            frame++;
        }
    }
   // printf("%f\n",tmp3/double(frame));
    rms[0] = sqrt(tmp3 / (double(frame)));
      rewind(wav);
      fclose (wav);
      
}

void loadWavWithDelay(char * filename, char * outFilename, const char *jsonFilename, char *treeDir, float gain, int frameAve, float thresh, int verbose,int delay,char* treeFilePrefix) 
{

    //    sprintf(str1, "trees/%s/levelClass", treeDir);
    //    sprintf(str2, "trees/%s/snrClass", treeDir);
     char str1[100] , str2[100];
    FILE * pFile;
     
            
    sprintf(str1,treeFilePrefix );
    DTree distTree;
    if (verbose ==1){printf("\nLoading Decision Trees");}
    distTree.readTextFilesTrees(str1);
    if (verbose ==1){printf("\nDone");}
    

    //http://www.sonicspot.com/guide/wavefiles.html
    //http://yannesposito.com/Scratch/en/blog/2010-10-14-Fun-with-wav/char 
    //char  filename[]="test16S.wav";
    //char  filename[]="iphone1.wav";
    float* maxL = (float*) malloc(sizeof (float)*1);
    maxL[0] = 1.5;
    float *rms = (float*) malloc(sizeof (float)*1);
    rms[0] = 0;
    wavRms(filename, maxL, rms);    
    double rmsD=(double)rms[0];

    FILE * wav;
    wav = fopen(filename, "r");
    //check openable
    int test = sizeof (header);
    if (fread(&header, sizeof (header), 1, wav) < 1) {
        fprintf(stderr, "Can't read input file header %s\n", filename);
        exit(1);
    }
    if (verbose ==1){
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
    }
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
     NoBlocks = datasize / (header.bits_per_sample / 8 * header.nochan);

    // allocate memory for time domain items
    double *windowD, *windowOverD, *windowPrevD;
    float *window;
    float *windowOver;
    float *windowPrev;
    float *last1s;
    float global_scale_max_over_1s = 0;
    float lastscale = 0;
    float scale =0;
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
    double rms1s = 0;
    double SpecFlux = 0;
    double SpecFlux2 = 0;

   // FILE * pFileRawData;
   // char ss[100];
   // sprintf(ss,"rawdata_delay_%i.txt",delay);     
//    remove("rawdata.txt");
//    pFileRawData = fopen("rawdata.txt", "w");
   // remove(ss);pFileRawData = fopen(ss, "w");
    
   // FILE *pFileMeanData;
    //remove("meandata.txt");
    //pFileMeanData = fopen("meandata.txt", "w");
   //    sprintf(ss,"meandata_delay_%i.txt",delay);  
  //  pFileMeanData = fopen(ss, "w");
    
   sprintf(str2,"%s_%i.txt",outFilename,delay);  
    //sprintf(str2, "%s.b", outFilename);
    remove(str2);
    pFile = fopen(str2, "w");


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
                     if (verbose ==1){
                    printf("Supported formats are 16 and 32 bit signed integer, and IEEE float little-endian - any number of channels");
                     }
                    break;
                }

                if (i>delay){
                double tmp3 = ((double) tmp2 / (double) header.nochan); //*  158489=  10^(104/20)

                window[wN] = last1s[N1s]; //last window before overwritten

                if (abs(tmp3) > global_scale_max_over_1s)
                    global_scale_max_over_1s = abs(tmp3);

                    last1s[N1s] = tmp3; // Overwrite

                // printf("%f \n",window[wN]);
             
                N1s++;
                wN++;
                if (N1s == (frameAve * WIN_N)) {
                    N1s = 0;                   
                    lastscale=global_scale_max_over_1s;
                    scale=global_scale_max_over_1s;
                    global_scale_max_over_1s = 0;
                }
                }
            }

            //every WIN_N samples do this
            if (wN == WIN_N ) {
                                

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
                
                rms1s = rms1s + scaleover_temp / ((double) frameAve ) + scale_temp / ((double) frameAve );

                scaleover_temp = sqrt(scaleover_temp);
                scaleprev_temp = sqrt(scaleprev_temp);
                scale_temp = sqrt(scale_temp);
                
                double maxNormOverD =0;
                double maxNormD=0;
                for (n = 0; n < WIN_N; n++) {

                    //                tdData=tdData./sqrt(mean((tdData.^2)));
                    //                scale=max(abs(x1s(seci,:)./sqrt(mean(x1s(seci,:).^2))))/rms_store(i);
                    windowOverD[n] = ((double) windowOver[n]) / scaleover_temp;
                    windowPrevD[n] = ((double) windowPrev[n]) / scaleprev_temp;
                    windowD[n] = ((double) window[n]) / scale_temp;
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
                //normed2rms
                //scale=max(normed2rms)/originalrms
                double scaleC=lastscale / scaleover_temp;
                gatewayFunctionDist(windowOverD,
                        WIN_N,
                        (double) header.fs,
                        scaleC,
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
            count_[0]=count_[0]/((double)WIN_N/(double)header.fs*2.0);
                for (n = 0; n < 255; n++) {
                    featstmp[n] = outMatrix[n];
                }
                //     'hist'
            //     'spectral_skewness_Mean'
            //     'spectral_kurtosis_Mean'
            //     'spectral_spectentropy_Mean'
            //     'spectral_roughness_Mean'
            //     'timbre_zerocross_Mean'
            //     'timbre_spectralflux_Mean'
                
                featstmp[255] = skewness[0];
                 featstmp[256] = kurt[0];
               featstmp[257] = ent[0];
                featstmp[258] = roughness[0];
                featstmp[259] = count_[0];
                featstmp[260] = (SpecFlux2);
        


            

               // for (n = 0; n < 261; n++) {
               //     fprintf(pFileRawData, "%f ", featstmp[n]);
               // }
               // fprintf(pFileRawData, "\n");



                for (n = 0; n < 255; n++)
                    outMatrix_[n] = outMatrix_[n] + outMatrix[n];
                centroid_[0] = centroid_[0] + centroid[0] / ((double) frameAve * 2.0);
                spread_[0] = spread_[0] + spread[0] / ((double) frameAve * 2.0);
                kurt_[0] = kurt_[0] + kurt[0] / ((double) frameAve * 2.0);
                ent_[0] = ent_[0] + ent[0] / ((double) frameAve * 2.0);
                roughness_[0] = roughness_[0] + roughness[0] / ((double) frameAve * 2.0);
                skewness_[0] = skewness_[0] + skewness[0] / ((double) frameAve * 2.0);
                count__[0] = count__[0]+count_[0]/ ((double) frameAve * 2.0);
                ;
                ;


                for (n = 0; n < WIN_N; n++) {
                    peaks[n] = 0;
                }
                //Overlapping window
                for (n = 0; n < 255; n++)
                    outMatrix[n] = 0;
                
                scaleC=lastscale / scale_temp;
                gatewayFunctionDist(windowD,
                        WIN_N,
                        (double) header.fs,
                        scaleC,
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
            count_[0]=count_[0]/((double)WIN_N/(double)header.fs*2.0);

                for (n = 0; n < 255; n++) {
                    featstmp[n] = outMatrix[n];
                }
                   featstmp[255] = skewness[0];
                 featstmp[256] = kurt[0];
               featstmp[257] = ent[0];
                featstmp[258] = roughness[0];
                featstmp[259] = count_[0];
                featstmp[260] = (SpecFlux2);

               // for (n = 0; n < 261; n++) {
               //     fprintf(pFileRawData, "%f ", featstmp[n]);
               // }
               // fprintf(pFileRawData, "\n");

                for (n = 0; n < 255; n++)
                    outMatrix_[n] = outMatrix_[n] + outMatrix[n];
                centroid_[0] = centroid_[0] + centroid[0] / ((double) frameAve * 2.0);
                spread_[0] = spread_[0] + spread[0] / ((double) frameAve * 2.0);
                kurt_[0] = kurt_[0] + kurt[0] / ((double) frameAve * 2.0);
                ent_[0] = ent_[0] + ent[0] / ((double) frameAve * 2.0);
                roughness_[0] = roughness_[0] + roughness[0] / ((double) frameAve * 2.0);
                skewness_[0] = skewness_[0] + skewness[0] / ((double) frameAve * 2.0);
                count__[0] = count__[0]+count_[0] / ((double) frameAve * 2.0);



                SpecFlux = 0;
                for (n = 0; n < WIN_N / 2 + 1; n++)
                    SpecFlux = SpecFlux + ((pow(spectrum[n] - lastOverspectrum[n], 2)));
                SpecFlux2 = SpecFlux2 + sqrt(SpecFlux) / ((double) frameAve * 2.0);
                //            printf("%f \n",sqrt(SpecFlux));
                SpecFlux = 0;
                for (n = 0; n < WIN_N / 2 + 1; n++)
                    SpecFlux = SpecFlux + ((pow(lastOverspectrum[n] - lastspectrum[n], 2)));
                SpecFlux2 = SpecFlux2 + sqrt(SpecFlux) / ((double) frameAve * 2.0);
                counter++;
                
                // Every )frameAve frames do, (not counting the 50% overlap)
                if ((N1s == 0) && (frame > frameAve)) {
                //if ((frame % (frameAve)) == 0 && (frame > frameAve)) {
                    double pmfmean = 0;
                    for (n = 0; n < 255; n++) {
                        pmfmean = pmfmean + outMatrix_[n];
                    }

                    for (n = 0; n < 255; n++) {
                        feats[n] = outMatrix_[n] / pmfmean;
                        //printf("%f \n",feats[n]);
                    }

                    feats[255] = skewness_[0];
                    feats[256] = kurt_[0];
                    feats[257] = ent_[0];
                    feats[258] = roughness_[0];
                    feats[259] = count__[0];
                    feats[260] = (SpecFlux2);

                    //                printf("%f \n",feats[260]);
                   // for (n = 0; n < 261; n++) {
                   //     fprintf(pFileMeanData, "%f ", feats[n]);
                    //    //printf("%f ",feats[n]);
                  //  }
                    int distLevel = distTree.decisionTreeFun(feats);
                    if(distLevel==0)
                        distLevel=5; //as 0 is a fault, probably an nan somewhere, assume full quality

                    //printf("frame centre %f Distortion Amount %i \n",(frame),5-distLevel);
                    fprintf(pFile, " %0.2f %0.2f %i \n", (double) (i-WIN_N*frameAve) / 44100.0, sqrt(rms1s)/rmsD, distLevel);
                    if (verbose ==1){
                    printf("Time (s) %0.2f RMS %0.2f Quality %0.0f \n", (double) (i  -WIN_N*frameAve) / 44100.0, sqrt(rms1s)/rmsD , ((double) distLevel-1)/4.0*100.0);
                    }
                    //fprintf(pFileMeanData, "\n");
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
                    rms1s = 0;

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
//    fclose(pFileRawData);
   // fclose(pFileMeanData);
    
}




int loadWav(char * filename, char * outFilename, const char *jsonFilename, char *treeDir, float gain, int frameAve, float thresh, int verbose,char* treeFilePrefix)
{
   
    char str1[100], str2[100] ,str3[100] , datastring[100000];
    FILE *pFileJSON;
    FILE * pFile;
   float* maxL = (float*) malloc(sizeof (float)*1);
    maxL[0] = 1.5;
    float *rms = (float*) malloc(sizeof (float)*1);
    rms[0] = 0;
    wavRms(filename, maxL, rms);    
    
    // NoBlocks = datasize / (header.bits_per_sample / 8 * header.nochan);

    double over = 0.25;
    double * timeStore = (double*) calloc ((int)(NoBlocks/(frameAve*WIN_N*over)),sizeof(double));
    double * rmsStore = (double*) calloc ((int)(NoBlocks/(frameAve*WIN_N*over)),sizeof(double));
    double * rmsStore2 = (double*) calloc ((int)(NoBlocks/(frameAve*WIN_N*over)),sizeof(double));
    double * QualStore = (double*) calloc ((double)(NoBlocks/(frameAve*WIN_N*over)),sizeof(double));
    double * QualStore2 = (double*) calloc ((double)(NoBlocks/(frameAve*WIN_N*over)),sizeof(double));
    int ii=0;
           for (ii=0;ii<(int)(NoBlocks/(frameAve*WIN_N*over));ii++){
               QualStore[ii]=1;
               QualStore2[ii]=1;
           }
       

    int  nn=0;
    int cc=0;
    ii=0;
    for (ii=0;ii<(int)(1/over);ii++){
        
        int delay=(int)((float)(WIN_N*frameAve)*(over*ii));
        loadWavWithDelay(filename, outFilename, jsonFilename, treeDir,  gain,  frameAve,  thresh,  verbose,delay ,treeFilePrefix);
        sprintf(str2,"%s_%i.txt",outFilename,delay);  
        pFile = fopen(str2, "r");
        nn=0;
        while (fgets(datastring, sizeof (datastring), pFile) != NULL) {
            cc++;
            float rms, t, Qual;
            sscanf(datastring, "%f %f %f", &t, &rms,
                &Qual);
            timeStore[nn]=(double)t;
            rmsStore[nn* (int)(1/over)+ii]=(double)rms;
           QualStore[nn* (int)(1/over)+ii]=(double)(Qual-1)/4;
           QualStore2[nn* (int)(1/over)+ii]=(double)(Qual-1)/4;

         // printf("%f %f %f \n", t, rms,Qual);
            nn++;
        }
    }
    
    sprintf(str2, "%s.b", outFilename);
    remove(str2);
    pFile = fopen(str2, "w");
    //Weighted  sum denom
    double weight=1;
    for (int jj=1;jj<(1/over);jj++){
       weight=weight+ (1-jj*over)*2;
    }
    for (ii=0;ii<(cc+(int)1/over);ii++)
    {
        if (ii>=((int)1/over))
        {            
            QualStore2[ii]=QualStore2[ii]/(weight);
            rmsStore2[ii]=rmsStore2[ii]/(weight);
            for (int jj=1;jj<(1/over);jj++){
            QualStore2[ii]=QualStore2[ii]+QualStore[ii-jj]*(1-jj*over)/(weight)+QualStore[ii+jj]*(1-jj*over)/(weight);
            rmsStore2[ii]=rmsStore2[ii]+rmsStore[ii-jj]*(1-jj*over)/(weight)+rmsStore[ii+jj]*(1-jj*over)/(weight);
            }


        } 
          if (ii<((int)1/over))  {
            QualStore2[ii]=QualStore[ii];
            rmsStore2[ii]=rmsStore[ii];
          }
        
           //printf("%f\n",QualStore2[ii]);
            double rms, t, Qual;
            t=(ii+1)*frameAve*WIN_N*over/(double)header.fs;
            Qual=(QualStore2[ii]*4);
            fprintf(pFile,"%f %f %f\n", t, rmsStore[ii],Qual);
            //printf("%f %f %f\n",t, rmsStore[ii],Qual);

    }
    fclose(pFile);

    // Global stats
    float averms = 0;
    float maxrms = 0;
    float minrms = 1000000;

    //float thresh = 0.2;
    float aveQual = 0;
    float maxQual = 5;
    float minQual = 0;

    float count_wF5 = 0;
    float count_wF4 = 0;
    float count_wF3 = 0;
    float count_wF2 = 0;
    float count_wF1 = 0;
    float count_wF0 = 0;


   int counter = 0;
    //  pFile2 = fopen ( 'tmp.txt' , "w");
    pFile = fopen(str2, "r");
    char mystring[100000];
    while (fgets(mystring, sizeof (mystring), pFile) != NULL) {

        float rms, t, Qual;
        sscanf(mystring, "%f %f %f", &t, &rms,
                &Qual);

        if (rms > maxrms)
            maxrms = rms;
        if (rms < minrms)
            minrms = rms;


        if (Qual > maxQual)
            maxQual = Qual;
        if (Qual < minQual)
            minQual = Qual;



        averms = averms + (rms);
        aveQual = aveQual  + Qual* rms;//; //%quality weighted by rms level

        Qual=4-Qual;
        if (Qual >= 0 && Qual <= 0.5)
            count_wF5++;
        if (Qual > 0.5 && Qual < 1.5)
            count_wF4++;
        if (Qual >= 1.5 && Qual < 2.5)
            count_wF3++;
        if (Qual >= 2.5 && Qual < 3.5)
            count_wF2++;
        if (Qual >= 3.5 && Qual < 4)
            count_wF1++;

   // printf( "%f, %0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f\n",Qual, count_wF1, count_wF2, count_wF3, count_wF4, count_wF5)    ;  
    counter++;
    }
    fclose(pFile);
    ;
    count_wF1 = count_wF1 / counter * 100.0;
    count_wF2 = count_wF2 / counter * 100.0;
    count_wF3 = count_wF3 / counter * 100.0;
    count_wF4 = count_wF4 / counter * 100.0;
    count_wF5 = count_wF5 / counter * 100.0;

    aveQual = (aveQual / 4) / averms;
    averms = averms / counter;

    // remove(outFilename);
    pFileJSON = fopen(jsonFilename, "w");
    //

    FILE * pFile2;
    pFile2 = fopen(outFilename, "w");
    pFile = fopen(str2, "r");
    fprintf(pFile2, "Distortion Detection - University of Salford - the Good Recording Project http://www.goodrecording.net \n\n", filename);
    fprintf(pFile2, "Distortion Analysis for input file %s\n\n", filename);
    fprintf(pFile2, "Distortion Statistics, %% number of frames with distortion detected at the following Quality Levels \n\n");
    fprintf(pFile2, "%% of time in each Degradation range.\n");
    fprintf(pFile2, "Bad,\tPoor,\tFair,\tGood,\tExcellent\n");
    fprintf(pFile2, "%0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f\n", count_wF1, count_wF2, count_wF3, count_wF4, count_wF5);
    fprintf(pFile2, "\n");
    fprintf(pFile2, "Distortion noise time history\n\n");
    fprintf(pFile2, "T(s)\t\tQuality Degradation(%%)\t RMS (normalised to global RMS) \n");
    //JSON global stats
    //printf("heeeerrrre %s\n",jsonFilename);
    fprintf(pFileJSON, "{\n \t\"Global Stats\":\n \t [ \n");
    fprintf(pFileJSON, "\t\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f,\t%0.1f\n\t\t],\n", count_wF1, count_wF2, count_wF3, count_wF4, count_wF5);
    fprintf(pFileJSON, "\t \"Time History\":\n\t [\n");
    int firstWin = 0;
    while (fgets(mystring, sizeof (mystring), pFile) != NULL) {

        float rms, t, qual;
        sscanf(mystring, "%f %f %f", &t, &rms,
                &qual);

        fprintf(pFile2, "%0.2f\t\t%0.2f\t\t\t\t\t\t%0.2f\n", t, 100 - (qual / 4)*100.0, rms);
        if (firstWin == 1)
            fprintf(pFileJSON, ",\n");

        fprintf(pFileJSON, "\t\t{\"Ts\": %0.2f, \"Te\": %0.2f, \"rms\": %0.2f, \"QDeg\": %0.2f}", t-frameAve*WIN_N*over/(double)header.fs,t, rms, 100 - (qual / 4)*100.0);
        firstWin = 1;
        counter++;
    }    
    fclose(pFile);
    pFile = fopen(str2, "r");
   ///////////////////////////// Average Stats
    fprintf(pFile2, "\nAverage Quality (RMS weighted ave)\n %0.1f\n", ((double)aveQual)*100.0);
    fprintf(pFileJSON, "\n\t],\n");
    fprintf(pFileJSON, "\t\"AveRmsWQual\":\n \t [ \n");
    fprintf(pFileJSON, "\t\t%0.1f\n\t\t],\n", ((double)aveQual)*100.0);

    
    ///////////////////////////// Find contiguous regions without distortion noise

    firstWin = 0;
    int lastclean = 0; // 
    int thisclean = 0; // 
    int initialised = 0;
    int first = 0;
    float win = 0;
    float start = 0;
    fprintf(pFile2, "\nDistortion free regions from - to (s) using a Threshold of %2.0f\n\n", thresh);
    fprintf(pFileJSON, "\t\"Distortion free regions\":\n \t [\n");
    float rmsw, tw, Qualw;
    while (fgets(mystring, sizeof (mystring), pFile) != NULL) {
        sscanf(mystring, "%f %f %f", &tw, &rmsw,
                &Qualw);
        tw = tw - frameAve*WIN_N*over/(double)header.fs;
        Qualw = (Qualw / 4)*100.0; //Quality
        thisclean = (int) (Qualw > thresh);

       // printf("%i %i %f %f\n", thisclean, lastclean, Qualw, thresh);

        if (thisclean == 1 && lastclean == 0) // previous frame was dist now its clean
        {
            start = tw;
        } else if (thisclean == 0 && lastclean == 1) // previous frame was clean now its dist
        {
            fprintf(pFile2, "%0.2f\t%0.2f\n", start, tw-frameAve*WIN_N*over/(double)header.fs);
            if (firstWin == 1)
                fprintf(pFileJSON, ",\n");
            fprintf(pFileJSON, "\t\t{\"s\": %0.2f, \"e\": %0.2f}", start, tw-frameAve*WIN_N*over/(double)header.fs);
            firstWin = 1;
        }
        lastclean = thisclean;
        counter++;
    }
    // finish up if last window is clean
    if (thisclean == 1 && lastclean == 1) // previous frame was clean now its dist
    {
        fprintf(pFile2, "%0.2f\t%0.2f\n", start, tw);
        if (firstWin == 1)
            fprintf(pFileJSON, ",\n");
        fprintf(pFileJSON, "\t\t{\"s\": %0.2f, \"e\": %0.2f}", start, tw+frameAve*WIN_N*over/(double)header.fs);
        firstWin = 1;
    }

    ///////////////////////////
    fclose(pFile);
    fclose(pFile2);
    fprintf(pFileJSON, "\n\t]\n}\n");
    fclose(pFileJSON);
    return 0;
    //remove(str2);

}



/*==========================================================
 * Dist features extraction
 *========================================================*/
/* $Revision: 1.1.10.3 $ */
#define M_PI 3.14159265358979323846
#include "math.h"
#include "distFeats.h"
#include "kiss_fft130/kiss_fft.h"

// The PMF
void histFunc(double *value, double *bin, int ncols,double scale)
{
    double max=0;
    double min=0;
    int i=0;
    double interval = (double)(2.0 ) / 254.0;
   
//   Reset histogram
    for ( i=0;i<255;i++)
        bin[i]=0;
    // histogram counting function uses int cast (which always rounds down so add 0.5)
    for ( i=0;i<ncols;i++)
    {
        int ii;
        ii=(int)(((value[i]/scale+1)/interval)+0.5);
        if (ii>=0 && ii<255)
            bin[ii]=bin[ii]+1;
//         printf ("%i,",(int) (((value[i]/scale+1)/interval)+0.5));
 //   bin[(int) ((value[i]/scale)/interval  +0.5) ]=bin[(int) ((value[i]/scale)/interval  +0.5) ]+1;
    }
 for ( i=0;i<255;i++)
        bin[i]= bin[i]/(double)ncols;        
//         bin[ (int)((value[i]/scale- (-1.0))/interval+0.5) ]=bin[ (int)((value[i]/scale- (-1.0))/interval+0.5) ]+1.0/(double)ncols;
    
}

//Zero crossings

void zeroCrossings(double *value, int crossings[], int ncols,double Fs)
{
    double last=0;  
    int i=0;    
    crossings[0]=0;

    for ( i=0;i<ncols;i++)
    {
      
        if(value[i]>0 && last<0)
            crossings[0]=crossings[0]+1;
        if(value[i]<0 && last>0)
            crossings[0]=crossings[0]+1;
        last=value[i];
    }
         //  crossings[0]=crossings[0]/((1024.0/Fs)/2);

    
    
}


// Compute the spectrum
void computeSpectrum(double *value, double *spectrum,int nfft)
{
    size_t buflen = sizeof(kiss_fft_cpx)*nfft;
    
    kiss_fft_cpx  * in = (kiss_fft_cpx*)malloc(buflen);
    kiss_fft_cpx  * out= (kiss_fft_cpx*)malloc(buflen);
    kiss_fft_cfg  cfg = kiss_fft_alloc(nfft,0,0,0);
    int k;
    double win;
    double min=10e8;
    for (k=0;k<nfft;++k) {
        win=0.54 - 0.46 * cos(2*M_PI*k/(nfft-1)); //hamming window
        in[k].r = value[k]*win;
        in[k].i = 0;
        
    }
    kiss_fft(cfg,in,out);
    
    for (k=0;k<(nfft/2+1);++k)
    {
        spectrum[k]= (double)(pow(pow(out[k].r,2) +  pow(out[k].i,2),0.5));  //scaling as MIR toolbox does
        if (spectrum[k]<min)
            min=spectrum[k];
    }
    
    for (k=0;k<(nfft/2+1);++k)
    {
        spectrum[k]=spectrum[k]+min;
    }
    free(in);
    free(out);
    free(cfg);
}

//  Spectral centroid
void computeSpecCentroid( double *spectrum,double* out, double * f,int nfft)
{
    //centrue=sum(f.*spec)/sum(spec)
    double spectrumSum = 0;int k;out[0]=0;
    for (k=0;k<(nfft/2+1);++k)
        spectrumSum=spectrum[k]+spectrumSum;
    for (k=0;k<(nfft/2+1);++k)
        out[0]=  out[0] + (f[k]*(spectrum[k]));
    out[0]=out[0]/spectrumSum;

}

void computSpecSkew( double *spectrum,double out[], double * f,int nfft,double *centroid, double *spread)
{
//     s = sum((f-centroid).^3.* (spectrum/sum(spectrum)) ) ./ spread.^3;
    double spectrumSum = 0;int k;out[0]=0;
    
    for (k=0;k<(nfft/2+1);++k)
        spectrumSum=spectrum[k]+spectrumSum;
    
    for (k=0;k<(nfft/2+1);++k)
            out[0]=  out[0] + (  pow(f[k] - centroid[0],3) *(spectrum[k]/spectrumSum))/pow(spread[0],3);
}
//  Spectral spread
void computeSpecSpread( double *spectrum,double out[], double * f, double *cen,int nfft)
{
    // sqrt( sum((f-cen(i)).^2.*spn/sum(spn))  );
    double spectrumSum = 0;int k;out[0]=0;
    for (k=0;k<(nfft/2+1);++k)
        spectrumSum=spectrum[k]+spectrumSum;
    for (k=0;k<(nfft/2+1);++k)
        out[0]=  out[0] + pow((f[k]-cen[0]),2) * spectrum[k]/spectrumSum;
    out[0]=pow(out[0],0.5);
}

//  Spectral kurt
void computeSpecKurt( double *spectrum,double out[], double * f, double *cen,double* spread,int nfft)
{
// kurttrue = sum((f-cen(i)).^4.*(spn/sum(spn))) ./ spread.^4;
    double spectrumSum = 0;int k; out[0]=0;
    for (k=0;k<(nfft/2+1);++k)
        spectrumSum=spectrum[k]+spectrumSum;
    for (k=0;k<(nfft/2+1);++k)
        out[0]=  out[0] + pow((f[k]-cen[0]),4) * (spectrum[k]/spectrumSum);
    out[0]=out[0] / pow(spread[0],4) ;
}

//  Spectral Entropy
void computeSpecEnt( double *spectrum,double out[],int nfft)
{
//  enttrue=- sum(spn.*log(spn + 1e-12))./log(size(spn,1));
    double min=10e8;
    double spectrumSum = 0; int k;    out[0]=0;
    for (k=0;k<(nfft/2+1);++k)
    {
        if(spectrum[k]<min)
            min=spectrum[k];
    }
    
    for (k=0;k<(nfft/2+1);++k)
        spectrumSum=(spectrum[k] - min/2)+spectrumSum;
    for (k=0;k<(nfft/2+1);++k)
        out[0]=  out[0] + ((spectrum[k]- min/2)/spectrumSum) * log((spectrum[k]- min/2)/spectrumSum+1e-12);
    out[0]=-out[0] / log(nfft/2+1) ;
}

// Find peaks in spectrum
void computeSpecPeaks( double *spectrum,double * out, double * f,int nfft)
{
    // Finds all peaks
    //
    double max;
    double * pks;
    double *    vlys;
    double cth,min;
    int k,j;
    double lvlnextyval,lvllastyval;
    
    pks = (double *)calloc((nfft/2+1),sizeof(double));
    vlys = (double *)calloc((nfft/2+1),sizeof(double));
    for (k=1;k<(nfft/2);++k)
    {
        if (spectrum[k]>spectrum[k-1] && spectrum[k]>spectrum[k+1])
            pks[k]=1;
    }
    
    // Find valleys in spectrum
    //
    max=0;
    min=10e8;
    for (k=1;k<(nfft/2);++k)
    {
        if (spectrum[k]<spectrum[k-1] && spectrum[k]<spectrum[k+1])
            vlys[k]=1;
        if (spectrum[k]>max)
            max=spectrum[k];
        if (spectrum[k]<min)
            min=spectrum[k];
        
    }
    cth= (max-min)*0.01;
    //Remove peaks out of range
    for (j=1;j<(nfft/2);++j)
    {
        //level of last valley
        lvllastyval=0;
        for (k=0;k<(j);++k)
        {
            if (vlys[k]>0)
                lvllastyval=spectrum[k];
        }
        
        //level of next valley
        lvlnextyval=0;
        for (k=j+1;k<(nfft/2+1);++k)
        {
            if (vlys[k]>0)
            {
                lvlnextyval=spectrum[k];
                break;
            }
        }
        
        if (pks[j]==1)
        {
            if ((spectrum[j]-lvllastyval)>=cth && (spectrum[j]-lvlnextyval)>=cth)
            {
                out[j]=1;
            }else
            {
                out[j]=0;
            }
            
        }
    }
    
    free(vlys);free(pks);
}


//  Roughness
void computeRoughness( double *spectrum,double * out,double *peaks, double * f,int nfft)
{
    int j,k,n;
    double f1,f2,b1,b2,xstar,s1,s2,s;    

    double min;
    double max;
   double * mulamp = (double *)calloc((nfft/2+1)*(nfft/2+1),sizeof(double));
    double *pd = (double *)calloc((nfft/2+1)*(nfft/2+1),sizeof(double));

 
    for (k=0;k<(nfft/2+1);++k)
    {
        if(spectrum[k]<min)
            min=spectrum[k];
        if(spectrum[k]>max)
            max=spectrum[k];   
//                 printf ("%f,", f[k]);

    }
    
    b1 = 3.51;
    b2 = 5.75;
    xstar = .24;
    s1 = .0207;
    s2 = 18.96;
    n=0;

    
    for (j=0;j<(nfft/2+1);j++)
    {
        for (k=0;k<(nfft/2+1);k++)
        {
            f1=f[j];f2=f[k];
            if (peaks[j]==1 && peaks[k]==1)
            {

                if(f1<f2)
                {
                    s = (xstar / (s1 * f1 + s2 ));
                }else if(f1>=f2)
                {
                    s = (xstar / (s1 * f2 + s2 ));
                }
                if (n>(nfft/2+1)*(nfft/2+1))
                    break;
                pd[n]=expf(-b1*s*fabs(f2-f1)) - expf(-b2*s*fabs(f2-f1));
//                 printf ("%f %f,", fabs(f2-f1),  expf(-b2*s*fabs(f2-f1))  );

                mulamp[n]=(spectrum[j])*(spectrum[k]); 
                n++;
            }           
        }
    }

    out[0]=0;
    for (j=0;j<n;j++)
    {
        out[0]=out[0]+(pd[j]*mulamp[j])/2.0;
    }
    free(mulamp);free(pd);
}



        
        
        /* The gateway function */
////void mexFunction( int nlhs, mxArray *plhs[],
 //       int nrhs, const mxArray *prhs[])
void gatewayFunctionDist(double *inMatrix,
        int ncols,
        double Fs,
        double scale,
    double *outMatrix, 
    double *spectrum,
    int* counter,
    double *centroid,                /* output matrix */
    double * spread,
    double *kurt,
    double *ent,
    double *peaks,
    double *roughness,
    double *skewness)
{
 
    /* 1xN input matrix */
    //size_t ncols;                   /* size of matrix */
    //double *Fst;
    //double *scalet;
    
 // double *outMatrix;              /* output matrix */
 //   double *spectrum;              /* output matrix */
 //   int* counter;                    /* output matrix */
 //   double *centroid;                /* output matrix */
 //   double * spread;
 //   double *kurt;
 //   double *ent;
 //   double *peaks;
 //   double *roughness;
 //   double *skewness;

    //Working vars
    //double scale;
    //double Fs;
    int k;
    double *f;
    /* create a pointer to the real data in the input matrix  */
   // inMatrix = mxGetPr(prhs[0]);
    //ncols = mxGetN(prhs[0]);
    //Fst= (double *)mxGetPr(prhs[1]);
    //scalet= (double *)mxGetPr(prhs[2]);

    //Fs=Fst[0];
   // scale=scalet[0];
    
    f = (double *)calloc(1024,sizeof(double));
    for (k=0;k<(ncols/2+1);++k)
    {
        f[k]=(((double)k)/((double)ncols/2))*((double)Fs/2);
    }
    
    /* call the computational routines */
    histFunc(inMatrix,outMatrix,ncols,scale);
//    printf("%i",counter[0]);
    zeroCrossings(inMatrix,counter,ncols,Fs);
    computeSpectrum(inMatrix,spectrum,ncols);
    computeSpecCentroid(spectrum,centroid, f,ncols);
    computeSpecSpread(spectrum,spread,f,centroid,ncols);
    computeSpecKurt(spectrum,kurt,f,centroid,spread,ncols);
    computeSpecEnt( spectrum,ent,ncols);
    computeSpecPeaks( spectrum,peaks, f,ncols);
    computeRoughness( spectrum,roughness,peaks, f,ncols);
    computSpecSkew( spectrum, skewness,  f, ncols, centroid,  spread);
//    printf("%f\n",counter[0]);
    free(f);
            
}

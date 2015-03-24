/* 
 * File:   disFeats.h
 * Author: ags056
 *
 * Created on 03 March 2015, 11:32
 */

#ifndef DISFEATS_H
#define	DISFEATS_H


#include <stdio.h>      /* Standard Library of Input and Output */
//#include <complex.h>    /* Standart Library of Complex Numbers */
 //void melSpec(int Fs,float * f, float * spec,int fftn,float *mfcc);
// void dctComp(float * data, int N,float *mfcc);
 

#ifdef	__cplusplus
extern "C" {
#endif

void histFunc(double *value, double *bin, int ncols,double scale);
void zeroCrossings(double *value, int crossings[], int ncols,double Fs);
void computeSpectrum(double *value, double *spectrum,int nfft);
void computeSpecCentroid( double *spectrum,double out[], double * f,int nfft);
void computSpecSkew( double *spectrum,double out[], double * f,int nfft,double *centroid, double *spread);
void computeSpecSpread( double *spectrum,double out[], double * f, double *cen,int nfft);
void computeSpecKurt( double *spectrum,double out[], double * f, double *cen,double* spread,int nfft);
void computeSpecEnt( double *spectrum,double out[],int nfft);
void computeSpecPeaks( double *spectrum,double * out, double * f,int nfft);
void computeRoughness( double *spectrum,double * out,double *peaks, double * f,int nfft);
void gatewayFunctionDist(double *inMatrix,
    int ncols,
    double Fs,
    double scale,
    double *outMatrix, 
    double *spectrum,
    int* counter,
    double *centroid,             
    double * spread,
    double *kurt,
    double *ent,
    double *peaks,
    double *roughness,
    double *skewness);

#ifdef	__cplusplus
}
#endif

#endif	/* DISFEATS_H */


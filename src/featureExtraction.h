/*
 * File:   featureExtraction.h
 * Author: ags056
 The MIT License (MIT)

Copyright (c) <2014> <Paul Kendrick>

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
 * Created on 10 April 2014, 16:19
 */

#ifndef FEATUREEXTRACTION_H
#define	FEATUREEXTRACTION_H

#include <stdio.h>      /* Standard Library of Input and Output */
//#include <complex.h>    /* Standart Library of Complex Numbers */
#include "aweighting.h"
 void melSpec(int Fs,float * f, float * spec,int fftn,float *mfcc);
 void dctComp(float * data, int N,float *mfcc);
 void aWeighting(float * f, float * a ,int Len, float * rms,int fftn );
 void aWeightingPreComp(float * f,int fftn);

#endif	/* FEATUREEXTRACTION_H */

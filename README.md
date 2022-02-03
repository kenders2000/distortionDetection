# Distortion and Clipping Detection in Audio Files
### University Of Salford Acoustics Research Centre

## Overview 

This program automatically analyses .wav files and detects regions where there may distortion due to overloading.  Distortion can degrade the quality of recordings.  

This tool could be used when there are a large number of recordings of a particular source, for example clips of recordings from mobile devices of an outdoor concert (user generated content).  This tool can quickly analyse the audio from each to determine samples with the least amount of distortion and therefore of the highest quality.

The distortion detection algorithm also provides time stamps identifying the regions which are free of distortion.  This can allow producers to quickly stitch together content from multiple sources, or flag regions where extra processing such as noise reduction is required.

Perceptual tests have shown that when distortion is present in speech music and other signals the HASQI (Hearing Aid Sound Quality index) is a good indicator of the quality of the recording for both normal and hearing impared listeners.  Therefore this detector works by trying to predict the HASQI which in turn indicates the possible degradation to the audio quality.

If you utilise this work in anyway please could cite our the following paper in your work:

Kendrick, Paul; Li, Francis; Fazenda, Bruno; Jackson, Iain; Cox, Trevor (2015). "Perceived Audio Quality of Sounds Degraded by Nonlinear Distortions and Single-Ended Assessment Using HASQI". Journal of the Audio Engineering Society 63 (9): 698–712.


## Usage

The program requires compiling using a C++ compiler and runs as a command line executable which analyses a wavfile and provides information as regards the level, location and effect on quality of any distortion that may be present.

Wav files must be 16 or 32 bit int, or 32 bit float PCM at 44.1 kHz. All channels are collapsed to one (average of all channels)

```
./distDet -i wav_filename -o output_filename  [-f 43] [-w 30] [-v 0]
     -i and -o are required parameters, they provide the input .wav filename and the output filename respectively. Specify the output as a text file, include path if required, the program will also output a json file with the same name in the same path.
     -f n, sets the number of frames used to produce the analysis window to n, windows are 1024 samples long.  The default is n = 43 frames, which is about 1 s
	 -w sets the threshold degradation in quality level for identifying noise free segments, default is 30(%), but this can be adjusted depending on the required quality level for the application
	 -v verbose 1 and diagnostic messages are printed, 0 and they are not
	 -j specifies the path and name of the json data file containing the json formated data.
```
   
### Output file description and example

The output text file will contain four sets of outputs: global statistics, distortion noise time history, and start and end points of clean regions. And the overall quality of the file, (weighted average of quality using the rms to weight) Using the following parameter set should provide the results presented.

```
./distDet "-v" "1" "-i" "test/02 BioD.wav" "-o" "test/02 BioD.txt"  "-f" "43" "-w" "25" "-j" "test/BioD.json"
```

#### Global statistics

This shows the percentage of frames which are either Bad, Poor, Fair, Good or Excellent quality.  These global statistics can be used to make a decision as regards the usability of the file.
Bad,	Poor,	Fair,	Good,	Excellent
12.5,	6.2,	0.0,	6.2,	75.0
The test output of the text file shows an example where 75.0% of frames are free of distortion noise, but 18.7% of the time the quality is fair or worse%.  (the json outputs the same data but formatted as to be readable by a json interpreter).

|Bad|Poor|Fair|Good|Excellent|
|0.0,|0.0,|9.7,|8.1,|82.3|

These statistics can be used to make quick judgements about the file regarding the distribution of distorted parts of the signal, but it is down to the required application to choose a suitable set of criteria for acceptability.

#### Distortion time history

The next set of parameters is the Distortion time history.  Three columns are presented: 

 1. The time at the end of the analysis window
 2. The predicted quality degradation which is computed from  the predicted HASQI
 3. The RMS level over the window normalised to the RMS level of the file  
 
```
|0.25		|0.00    					|0.94|
|0.50		|0.00	    					|0.75|
|0.75		|0.00						|0.77|
|1.00		|0.00						|0.74|
|1.25		|3.12						|0.85|
|1.50		|4.69						|0.90|
|1.75		|6.25						|0.92|
|2.00		|4.69						|0.93|
|2.25		|3.12						|0.93|
|2.50		|1.56						|0.91|
|2.75		|0.00						|0.89|
|3.00		|0.00						|0.87|
|3.24		|0.00						|0.83|
|3.49		|0.00						|0.78|
|3.74		|0.00						|0.71|
|3.99		|0.00						|0.62|
|4.24		|0.00						|0.53|
|4.49		|0.00						|0.53|
|4.74		|0.00						|0.66|
|4.99		|0.00						|0.98|
|5.24		|4.69						|1.48|
|5.49		|15.62						|2.02|
|5.74		|26.56						|2.52|
|5.99		|37.50						|2.87|
|6.24		|39.06						|3.05|
|6.49		|34.38						|3.12|
|6.74		|35.94						|3.13|
|6.99		|39.06						|3.08|
|7.24		|46.88						|2.98|
|7.49		|48.44						|2.79|
|7.74		|37.50						|2.46|
|7.99		|25.00						|2.06|
|8.24		|12.50						|1.65|
|8.49		|7.81						|1.31|
```
Etc.. Table is longer.. 
The average quality is the average predicted quality level (HASQI) but each quality value for each window is weighted by the rms level.

#### Average Quality (RMS weighted ave)
```
|82.9|
 ```
 
#### Start and end points distortion free sections -

Finally a set of times is returned which represent the start and end points, in seconds, of regions in the recording which are free of distortion.  So the following example does not show are regions above 75% quality degredation.


An intermediate file is also created which contains the raw classes for sanity checking.

## The Good Recording Project - overview

The University of Salford Acoustics Research Centre is carrying out a project entitled [The Good Recording Project](http://www.goodrecording.net/), investigating how recording errors can affect the perceived quality of audio.  In particular with regards to content recorded on mobile device such as mobile phones, video cameras and other such media often known as "user generated content".  This wind noise detector was developed as part of an [EPSRC funded project](http://gow.epsrc.ac.uk/NGBOViewGrant.aspx?GrantRef=EP/J013013/1).

It is licensed under the MIT license and uses the [KISS_fft library](http://sourceforge.net/projects.kissfft/) which is licensed under the revised BSD license. You are therefore free to use this in pretty much any project.  If you do use this work in any project [an acknowledgement would be appreciated](http://www.salford.ac.uk/computing-science-engineering/subjects/acoustics-audio-and-video).

The software was developed at the University of Salford and uses a bagged decision tree to analyse MFCCs to discover if an audio frame contains microphone wind noise, it is part of the Good Recording Project at the University of Salford Acoustic Research Centre. 


## How does it work?


For any comments or questions, please email kenders2000@gmail.com


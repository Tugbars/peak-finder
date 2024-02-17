#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include <stdbool.h>
#include "mqs_def.h"

struct peak {
    float magnitude;
    int location;
    float prominence;
};

#ifdef __cplusplus
extern "C"
{
#endif

    //extern struct peak mes_peakFinder(const int ndat, MqsRawDataPoint_t * input, bool includeEndpoints, float extrema); //mes_sweep.c'de callanmasi lazim yalniz bunun

    //extern struct peak calc_peakFinder(const int ndat, MqsRawDataPoint_t* input, MqsRawDataPoint_t* input2, bool includeEndpoints, float extrema);

    //extern struct peak calc_OldpeakFinder(const int ndat, MqsRawDataPoint_t* input, MqsRawDataPoint_t* input2, bool includeEndpoints, float extrema);

#ifdef __cplusplus
}
#endif
const float EPS = 2.2204e-16f;

/*
    Inputs
    x0: input signal
    extrema: 1 if maxima are desired, -1 if minima are desired
    includeEndpoints - If true the endpoints will be included as possible extrema otherwise they will not be included
    Output
    peakInds: Indices of peaks in x0
*/
//include endpoints feature has been deactivated.

   peak findPeaks(std::vector<float> x0, std::vector<peak>& peakHolder, std::vector<peak>& minpeakHolder, bool includeEndpoints = true, float extrema = 1.0);


#endif
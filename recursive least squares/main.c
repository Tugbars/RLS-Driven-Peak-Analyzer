#include "peakAnalysis.h"
#include "sliding_window_analysis.h"
#include "rls_polynomial_regression.h"
#include "trend_detection.h"
#include "statistics.h"
#include "mqs_def.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>


void myCallbackFunction(void) {
    if (boundaryErrorOccurred) {
        printf("Callback executed: State machine interrupted due to boundary error.\n");
        boundaryErrorOccurred = false;  // Reset the persistent flag after handling the error
    } else {
        printf("Callback executed: State machine returned to SWP_WAITING.\n");
    }
}

void PrepareBaseSweep(MesSweep_t *const sweep, MqsRawDataSet_t *const data)
{
  MesSweepApplySetupDefaults(&data->base);   // sets the parameters up such as increment, starting frequency etc.
  data->base.startFrequency = 210;
  MesSweepSetupAndClear(sweep, &data->base, data->data); // memset
}

int main() {
    
    PrepareBaseSweep(&rawBaseSweep, rawData); 
    currentRawSweep = &rawBaseSweep;
    
    int start_index = 162;
    startSlidingWindowAnalysis(currentRawSweep, start_index, myCallbackFunction); 
   
    return 0;
}


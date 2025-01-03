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

/**
 * @brief Prints the Savitzky-Golay window interval details.
 *
 * This function calls `get_savgol_window_interval`, retrieves the interval details,
 * and prints them in a formatted manner for debugging or logging purposes.
 */
void print_savgol_window_interval(void) {
    SavgolWindowInterval interval = get_savgol_window_interval();
    
    if (interval.analysisIntervalDistance == -1) {
        printf("Invalid analysis interval. Start or end index not set.\n");
        return;
    }
    
    printf("Savitzky-Golay Window Interval:\n");
    printf("  Adjusted Buffer Index: %d\n", interval.adjustedBufferIndex);
    printf("  Interval Distance: %d\n", interval.analysisIntervalDistance);
    printf("--------------------------------------\n");
}

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
    
    int start_index = 152; // 270 hata yapıyor
    startSlidingWindowAnalysis(currentRawSweep, start_index, myCallbackFunction); 
    
    //print_savgol_window_interval();
   
    return 0;
}

//too high negative value, do something about it.

/**
 * @file sliding_window_analysis.h
 * @brief Header file for sliding window analysis of phase angle data.
 *
 * This header defines the structures, enums, and function prototypes used for performing
 * sliding window analysis on phase angle data collected from an impedance analyzer.
 * The analysis is used to detect trends, peaks, and to adjust the analysis window dynamically.
 *
 */

#ifndef SLIDING_WINDOW_ANALYSIS_H
#define SLIDING_WINDOW_ANALYSIS_H

#include <stdint.h>
#include <stdbool.h>

// Dependencies
#include "peakAnalysis.h"
#include "buffer_manager.h"
#include "trend_detection.h"
#include "mes_buffers.h"

/**
 * @typedef Callback_t
 * @brief Typedef for a callback function used in sliding window analysis.
 *
 * This callback is invoked upon completion or at specific stages of the sliding window analysis.
 */
typedef void (*Callback_t)(void);

extern bool boundaryErrorOccurred;

extern bool undecided_error_flag;

/**
 * @enum SwpState_t
 * @brief Enumeration of states in the sliding window analysis state machine.
 *
 * This enum defines the various states that the sliding window analysis can be in.
 */
typedef enum {
    SWP_INITIAL_ANALYSIS,
    SWP_SEGMENT_ANALYSIS,
    SWP_UPDATE_BUFFER_DIRECTION,
    SWP_UNDECIDED_TREND_CASE,
    SWP_PEAK_CENTERING,
    SWP_PEAK_FINDING_ANALYSIS,
    SWP_PEAK_TRUNCATION_HANDLING,  
    SWP_EXPAND_ANALYSIS_WINDOW, 
    SWP_WAITING,
    SWP_STATE_LAST  // This should always be last
} SwpState_t;

/**
 * @struct SlidingWindowAnalysisContext
 * @brief Context structure for the sliding window analysis.
 *
 * This structure holds the current state and data required for performing the sliding window analysis.
 */
typedef struct {
    PeakPosition direction;       /**< Current direction based on peak analysis (LEFT, RIGHT, ON_PEAK, UNDECIDED). */
    Callback_t callback;          /**< Callback function to be invoked during the analysis. */
    bool isTruncatedLeft;         /**< Indicates if the left side of the sliding window is truncated. */
    bool isTruncatedRight;        /**< Indicates if the right side of the sliding window is truncated. */
    uint16_t peakIndex;      /**< Flag to indicate if the peak is near the boundary of the sliding window. */
    const double* phaseAngles;    /**< Pointer to the array of phase angles to be analyzed. */
    uint16_t phase_angle_size;    /**< Size of the phase angles array. */
} SlidingWindowAnalysisContext;

/**
 * @struct SavgolWindowInterval
 * @brief Struct to hold the adjusted buffer index and analysis interval distance.
 */
typedef struct {
    int adjustedBufferIndex; /**< The adjusted buffer index for the analysis start index. */
    int analysisIntervalDistance; /**< The absolute distance between analysis start and end indices. */
} SavgolWindowInterval;

/**
 * @brief Starts the sliding window analysis on the provided phase angle data.
 *
 * This function initializes and begins the sliding window analysis using the provided phase angle data.
 * It sets up the analysis context and invokes the analysis state machine, which progresses through
 * various states to detect trends, peaks, and adjust the analysis window accordingly.
 *
 * @param phaseAngles       Pointer to the array of phase angle data.
 * @param phase_angle_size  The size (number of elements) of the phase angles array.
 * @param callback          Callback function to be called during the analysis (optional).
 *
 * @note
 * - The phaseAngles array should contain the phase angle measurements collected from the impedance analyzer.
 * - The callback function can be used to perform actions upon certain events or completion of the analysis.
 */
void startSlidingWindowAnalysis(MesSweep_t *sweep, int start_index, Callback_t callback);

/**
 * @brief Sets the boundary error flag in the currentStatus.
 *
 * This function allows external modules (like buffer_manager) to set the boundary error flag in the sliding window
 * analysis state machine, which will trigger a state transition to SWP_WAITING.
 *
 * @param flag The value to set for the boundary error flag (1 for error, 0 to clear).
 */
void set_boundary_error_flag(uint8_t flag);


/**
 * @brief Computes and returns the analysis interval information for Savitzky-Golay filtering.
 *
 * This function combines the adjusted buffer index and the absolute distance
 * between the analysis start and end indices into a struct.
 *
 * @return SavgolWindowInterval A struct containing the adjusted buffer index and
 *         the analysis interval distance. If indices are invalid, returns -1 for both fields.
 */
SavgolWindowInterval get_savgol_window_interval(void);


#endif // SLIDING_WINDOW_ANALYSIS_H
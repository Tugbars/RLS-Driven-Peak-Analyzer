/**
 * @file peakAnalysis.c
 * @brief Implements peak significance analysis functions using MAD-based statistics.
 *
 * This file provides functions to analyze the significance of detected peaks within gradient trends.
 * It utilizes robust statistical measures to ensure accurate and reliable peak detection.
 */

#include "peakAnalysis.h"
#include "rls_polynomial_regression.h"  // Ensure RLS functions are accessible
#include "trend_detection.h"            // Include trend detection structures and functions
#include "statistics.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <stdio.h>
#include <stdarg.h>  // For va_list, va_start, va_end

// Define minimum trend length for analysis
#define MIN_TREND_LENGTH 5  // Adjust as necessary

/** ide
 * @def SIGNIFICANT_PEAK_THRESHOLD
 * @brief Threshold for maximum gradient to consider a peak as significant.
 */
#define SIGNIFICANT_PEAK_THRESHOLD 4.0

/** 
 * @def PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT
 * @brief Minimum number of consecutive points required to consider a trend as consistent during peak verification.
 */
#define PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT 5

/** 
 * @def PEAK_VERIFICATION_ALLOWABLE_INCONSISTENCY_COUNT
 * @brief Maximum number of consecutive minor fluctuations allowed during trend verification before considering the trend as broken.
 */
#define PEAK_VERIFICATION_ALLOWABLE_INCONSISTENCY_COUNT 1

/** 
 * @def SLOPE_THRESHOLD
 * @brief Threshold for second-order gradient to validate a peak.
 */
#define SLOPE_THRESHOLD 2.0

/** 
 * @def START_INDEX_OFFSET
 * @brief Number of positions to subtract from the start index for gradient tracking.
 */
#define START_INDEX_OFFSET 3

#define CENTERING_RATIO 2;  /**< Adjust this ratio to control centering aggressiveness */

#define NEGATIVE_GRADIENT_THRESHOLD (-1.0)

// Define debug levels
#define DEBUG_LEVEL 0  // Set to 0, 1, 2, or 3 to enable different levels of debugging

#if DEBUG_LEVEL >= 1
    #define DEBUG_PRINT_1(fmt, ...) debug_print(1, fmt, ##__VA_ARGS__)
#else
    #define DEBUG_PRINT_1(fmt, ...) do {} while(0)
#endif

#if DEBUG_LEVEL >= 2
    #define DEBUG_PRINT_2(fmt, ...) debug_print(2, fmt, ##__VA_ARGS__)
#else
    #define DEBUG_PRINT_2(fmt, ...) do {} while(0)
#endif

#if DEBUG_LEVEL >= 3
    #define DEBUG_PRINT_3(fmt, ...) debug_print(3, fmt, ##__VA_ARGS__)
#else
    #define DEBUG_PRINT_3(fmt, ...) do {} while(0)
#endif

#if DEBUG_LEVEL > 0
/**
 * @brief Prints debug messages based on the debug level.
 *
 * This function handles the actual printing of debug messages. It's separated from the macros to allow
 * for more complex debug handling in the future, such as logging to files or other outputs.
 *
 * @param level The debug level of the message.
 * @param fmt The format string.
 * @param ... The variable arguments corresponding to the format string.
 */
static void debug_print(int level, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    // Optionally, prepend debug messages with level information
    // printf("[DEBUG LEVEL %d] ", level);
    vprintf(fmt, args);
    va_end(args);
}
#endif

/**
 * @brief Analyzes gradient trends to determine peak characteristics.
 *
 * This function evaluates the gradient trends detected in a dataset to determine:
 * - Whether the current data window is centered on a peak (`isOnPeak`).
 * - Whether there is a significant peak within the data (`isSignificantPeak`).
 * - The direction to move the data window based on the trends (`moveDirection`).
 *
 * The analysis is based on the following criteria:
 * 
 * **1. Check if On a Peak:**
 *    - Calculates the lengths of both the increasing and decreasing trends.
 *    - Compares these lengths against a minimum trend length threshold (`MIN_TREND_LENGTH`).
 *    - If both trends have lengths greater than or equal to `MIN_TREND_LENGTH`, it sets `isOnPeak` to `true`, indicating that the window is centered on a peak.
 *
 * **2. Check for a Significant Peak:**
 *    - Evaluates the maximum gradient value (`maxValue`) from the increasing trend.
 *    - Compares `maxValue` against a significant peak threshold (`SIGNIFICANT_PEAK_THRESHOLD`).
 *    - If `maxValue` is greater than or equal to `SIGNIFICANT_PEAK_THRESHOLD`, it sets `isSignificantPeak` to `true`, indicating the presence of a significant peak.
 *
 * **3. Determine if Centering is Needed:**
 *    - If a significant peak is detected (`isSignificantPeak` is `true`) but the data window is not centered on the peak (`isOnPeak` is `false`)
 *    - This suggests that the window should be adjusted to better align with the peak.
 *
 * **4. Set Move Direction:**
 *    - Retrieves the move direction from the provided `trendResult`.
 *    - The move direction indicates how the data window should shift based on the detected trends:
 *        - `MOVE_RIGHT`: Shift the window to the right towards increasing trends.
 *        - `MOVE_LEFT`: Shift the window to the left towards decreasing trends.
 *        - `ON_THE_PEAK`: The window is centered on the peak; no movement needed.
 *        - `NO_TREND`: No dominant trend detected; movement is undecided.
 *
 * **Constants Used:**
 * - `MIN_TREND_LENGTH`: Defines the minimum length for a trend to be considered significant. Avoids using magic numbers.
 * - `SIGNIFICANT_PEAK_THRESHOLD`: Defines the threshold for the `maxValue` to consider a peak as significant.
 *
 * @param trendResult Pointer to the `GradientTrendResultAbsolute` structure containing the trend detection results.
 *
 * @return `PeakAnalysisResult` structure containing:
 * - `isOnPeak`: Indicates if the data window is centered on a peak.
 * - `isSignificantPeak`: Indicates if a significant peak is detected.
 * - `moveDirection`: The direction to move the data window based on trend analysis.
 */
PeakAnalysisResult analyzeGradientTrends(const GradientTrendResultAbsolute *trendResult) { // ON PEAK BULMA YONTEMI DETERMINE MOVE DIRECTIONDA COK KOTU. OLMAMASI DAHA IYI OLABILIR. 
    DEBUG_PRINT_1("Starting analysis of gradient trends.\n");

    // Initialize PeakAnalysisResult
    PeakAnalysisResult result = {
        .isOnPeak = false,
        .isSignificantPeak = false,
        .isWindowCentered = false,
        .moveDirection = UNDECIDED  // Default value from PeakPosition enum
    };

    // Initialize variables to hold trend lengths
    uint16_t increaseTrendLength = 0;
    uint16_t decreaseTrendLength = 0;
    result.moveDirection = trendResult->moveDirection; 

    // Calculate the length of the increasing trend
    if (trendResult->absoluteIncrease.startIndex != UINT16_MAX && trendResult->absoluteIncrease.endIndex != UINT16_MAX) {
        increaseTrendLength = trendResult->absoluteIncrease.endIndex - trendResult->absoluteIncrease.startIndex + 1;
        DEBUG_PRINT_2("Increase trend length: %u\n", increaseTrendLength);
    } else {
        DEBUG_PRINT_2("No valid increase trend detected.\n");
    }

    // Calculate the length of the decreasing trend
    if (trendResult->absoluteDecrease.startIndex != UINT16_MAX && trendResult->absoluteDecrease.endIndex != UINT16_MAX) {
        decreaseTrendLength = trendResult->absoluteDecrease.endIndex - trendResult->absoluteDecrease.startIndex + 1;
        DEBUG_PRINT_2("Decrease trend length: %u\n", decreaseTrendLength);
    } else {
        DEBUG_PRINT_2("No valid decrease trend detected.\n");
    }

    // Check for a significant peak first
    double maxValue = trendResult->absoluteIncrease.maxValue;  // Maximum gradient in the increase trend
    DEBUG_PRINT_2("Maximum value in increase trend: %.6f\n", maxValue);
    if (maxValue >= SIGNIFICANT_PEAK_THRESHOLD) {
        result.isSignificantPeak = true;
        result.moveDirection = ON_PEAK;
        printf("-->A significant peak is detected.\n");

        // If a significant peak is detected, check if the window is centered on the peak
        if (increaseTrendLength >= MIN_TREND_LENGTH && decreaseTrendLength >= MIN_TREND_LENGTH) {
            result.isOnPeak = true;
            printf("-->Data window is centered on the significant peak.\n");
        } else {
            DEBUG_PRINT_1("Significant peak detected but data window is not centered on it.\n");
        }
    } else {
        DEBUG_PRINT_1("No significant peak detected. Therefore, the window cannot be on a peak.\n");
    }

    // Determine if centering is needed
    if (result.isSignificantPeak && !result.isOnPeak) {
        result.isWindowCentered = false;
        printf("-->Centering is needed to align with the peak.\n");
    } else {
        DEBUG_PRINT_2("No centering needed.\n");
    }

    // Transfer moveDirection from trendResult to result, else ile girilmeli. 
                                                           // UNDECIDED IS RECEIVED HERE EVEN THOUGH PEAK VERIFIERS ARE RETURNING TRUE.
    printf("Move direction set to %s.\n",
                 (result.moveDirection == RIGHT_SIDE) ? "RIGHT_SIDE" :
                 (result.moveDirection == LEFT_SIDE) ? "LEFT_SIDE" :
                 (result.moveDirection == ON_PEAK) ? "ON_PEAK" :
                 (result.moveDirection == NEGATIVE_UNDECIDED) ? "NEGATIVE_UNDECIDED" : "UNDECIDED");

    DEBUG_PRINT_1("Completed analysis of gradient trends.\n");
    return result;
}


/**
 * @brief Performs peak analysis by identifying gradient trends and analyzing them.
 *
 * This function wraps the process of calling `identifyGradientTrends` to obtain the trend results
 * and then `analyzeGradientTrends` to analyze the trends. It simplifies the usage by encapsulating
 * both steps and returning a `PeakAnalysisResult`.
 *
 * **Additional Functionality:**
 *
 * If a peak or significant increase is detected (`isOnPeak` or `isSignificantPeak` is `true`),
 * the function performs an additional check:
 *
 * - Uses the interval from the significant increase (`trendResult.absoluteIncrease.startIndex` to `endIndex`).
 * - Adjusts the start index to 3 positions before `startIndex`, ensuring it doesn't go below zero.
 * - Calls `trackGradients` over this adjusted interval to calculate second-order gradients.
 * - Checks if any of the second-order gradients exceed a threshold (e.g., -1.0 or -2.0).
 * - If they do, it sets `isValidPeak` to `true`, indicating it's a valid peak based on curvature.
 *
 * If neither `isOnPeak` nor `isSignificantPeak` is `true`, it calls `identifySplitMoveDirection`
 * to verify the `moveDirection`.
 *
 * @param data             Array of data points (MqsRawDataPoint_t) to analyze.
 * @param dataLength       The total length of the data array.
 * @param startIndex       The starting index in the data array.
 * @param analysisLength   The number of data points to include in the analysis window.
 * @param degree           The degree of the polynomial used for regression.
 * @param gradientOrder    The order of the gradient (first-order or second-order).
 *
 * @return PeakAnalysisResult The result of the peak analysis.
 */
PeakAnalysisResult performPeakAnalysis(
    const MqsRawDataPoint_t *data,
    uint16_t dataLength,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree,
    GradientOrder gradientOrder
) {
    DEBUG_PRINT_1("Starting peak analysis.\n");

    // Initialize structures to hold results
    GradientTrendResultAbsolute trendResult;
    GradientCalculationResult   gradCalcResultInc;
    GradientCalculationResult   gradCalcResultDec;

    initializeGradientTrendResultAbsolute(&trendResult);
    memset(&gradCalcResultInc, 0, sizeof(GradientCalculationResult));
    memset(&gradCalcResultDec, 0, sizeof(GradientCalculationResult));

    DEBUG_PRINT_2("GradientTrendResultAbsolute initialized.\n");
    DEBUG_PRINT_2("Calling identifyGradientTrends with startIndex=%u, analysisLength=%u, degree=%u, gradientOrder=%d\n",
           startIndex, analysisLength, degree, gradientOrder);

    // Step 1: Perform Trend Detection
    identifyGradientTrends(
        data,
        startIndex,
        analysisLength,
        degree,
        gradientOrder,
        TREND_TYPE_GRADIENT_TRENDS,
        &trendResult,
        &gradCalcResultInc,
        &gradCalcResultDec
    );

    DEBUG_PRINT_3("identifyGradientTrends completed.\n");

    // Step 2: Analyze Gradient Trends
    PeakAnalysisResult peakAnalysisResult = analyzeGradientTrends(&trendResult);
    DEBUG_PRINT_2("Gradient trends analyzed.\n");

    // Initialize extra fields
    peakAnalysisResult.isValidPeak       = false;  // how is this defined?
    peakAnalysisResult.isWindowCentered  = false;

    // Step 3: If a peak or significant increase is detected, do additional checks
    if (peakAnalysisResult.isOnPeak || peakAnalysisResult.isSignificantPeak) {
        DEBUG_PRINT_1("Peak or significant increase detected. Proceeding with additional checks.\n");

        // We'll compute second-order gradients over the entire analysis window
        GradientCalculationResult secondOrderGradients;
        memset(&secondOrderGradients, 0, sizeof(GradientCalculationResult));

        uint16_t gradientStartIndex = startIndex;
        uint16_t gradientLength     = analysisLength;

        DEBUG_PRINT_3("Calculating second-order gradients over full window: startIndex=%u, length=%u\n",
                      gradientStartIndex, gradientLength);

        // Ensure we have enough data
        if (gradientStartIndex + gradientLength <= dataLength) { // not necessary
            // Step 3a: Calculate second-order gradients
            trackGradients(
                data,                     // MqsRawDataPoint_t array
                gradientLength,
                gradientStartIndex,
                degree,
                calculate_second_order_gradient,
                &secondOrderGradients
            );

            DEBUG_PRINT_3("Second-order gradients calculated.\n");

            // Step 3b: Inspect the last ~5 values
            int checkCount = 0;
            int checkStart = (secondOrderGradients.size >= 5)
                           ? (secondOrderGradients.size - 5)
                           : 0;
            DEBUG_PRINT_3("Checking the last 5 second-order gradients for threshold.\n");

            for (int i = checkStart; i < (int)secondOrderGradients.size; ++i) {
                double g = secondOrderGradients.gradients[i];
                DEBUG_PRINT_3("Gradient[%d] = %.6f\n", i, g);
                if (g < -1.0) { // or -2.0, or whichever threshold you like
                    checkCount++;
                    DEBUG_PRINT_3("Gradient[%d] is below threshold.\n", i);
                }
            }

            if (checkCount > 0) {
                // At least one gradient is below the threshold => window is sufficiently "centered"
                peakAnalysisResult.isWindowCentered = true;
                DEBUG_PRINT_1("Window is centered enough based on second-order gradients.\n");

                // Step 3c: Verify peak validity across entire window
                DEBUG_PRINT_3("Verifying peak validity using second-order gradient.\n");
                QuadraticPeakAnalysisResult verificationResult = verifyPeakValidity(
                    data, // also updated to accept MqsRawDataPoint_t if needed
                    dataLength,
                    startIndex,
                    startIndex + analysisLength - 1,
                    degree
                );

                peakAnalysisResult.isValidPeak = verificationResult.peak_found;
                if (verificationResult.peak_found) {
                    DEBUG_PRINT_1("Peak at index %u is valid (2nd order gradient).\n",
                                  verificationResult.peak_index);
                } else {
                    DEBUG_PRINT_1("Peak at index %u failed validity check.\n",
                                  verificationResult.peak_index);
                }

                // Optionally handle truncation flags
                if (verificationResult.is_truncated_left) {
                    DEBUG_PRINT_2("Peak at %u has left truncation.\n", verificationResult.peak_index);
                }
                if (verificationResult.is_truncated_right) {
                    DEBUG_PRINT_2("Peak at %u has right truncation.\n", verificationResult.peak_index);
                }

                // If window is centered and peak is valid, no centering needed
                peakAnalysisResult.isWindowCentered = true;
            } else {
                // Not centered enough => might need further centering
                peakAnalysisResult.isWindowCentered  = false;
                DEBUG_PRINT_1("Window is not centered enough; centering may be needed.\n");
            }
        } else {
            DEBUG_PRINT_1("Insufficient data to compute second-order gradients over the full window.\n");
            peakAnalysisResult.isWindowCentered = false;  // not necessary. 
        }
    }

    DEBUG_PRINT_1("Peak analysis completed.\n");
    return peakAnalysisResult;
}


/**
 * @brief Verifies the detected peak by checking for consistent trends on both sides.
 *
 * This function verifies whether a detected peak is a true peak by analyzing the second-order gradients
 * calculated using polynomial regression. Specifically, it checks for a consistent increasing trend on the left
 * side and a consistent decreasing trend on the right side of the peak.
 *
 * @param second_order_gradients Array containing the precomputed second-order gradients.
 * @param peak_index The index of the detected peak within the second_order_gradients array.
 * @param start_index The starting index in the original data array corresponding to the first element of second_order_gradients.
 * @return QuadraticPeakAnalysisResult Struct containing the verification result and truncation flags.
 */
/*************************************************************************
 * @brief Verifies the detected peak by checking for consistent trends on both sides.
 *
 * This function:
 *   1) Checks the left side for consistent positive second-order gradients,
 *      allowing a limited number of inconsistencies.
 *   2) Checks the right side for consistent negative second-order gradients,
 *      also allowing limited inconsistencies.
 *   3) Flags truncation if we run out of data on either side.
 *   4) Declares the peak verified only if both sides pass the criteria.
 *************************************************************************/
QuadraticPeakAnalysisResult verify_quadratic_peak(
    const double *second_order_gradients,
    uint16_t length,
    uint16_t peak_index,
    uint16_t start_index
) {
    QuadraticPeakAnalysisResult result = {
        .peak_found = false,
        .peak_index = (start_index + peak_index),
        .is_truncated_left = false,
        .is_truncated_right = false
    };

    /*************************************************************************
     * Variable Declarations
     *************************************************************************/
    uint16_t left_trend_count  = 0;
    uint16_t right_trend_count = 0;
    uint16_t inconsistency_count_left  = 0;
    uint16_t inconsistency_count_right = 0;

    /*************************************************************************
     * Debug
     *************************************************************************/
    DEBUG_PRINT_1("Verifying peak at absolute index %u.\n", result.peak_index);

    /*************************************************************************
     * Left-Side Check: look for consistent POSITIVE gradients.
     * We'll move left from `peak_index` down to 0 or until we exceed 
     * allowable inconsistencies.
     *************************************************************************/
    for (int i = peak_index; i > 0; --i) {
        double grad = second_order_gradients[i - 1];
        if (grad > 0) {
            left_trend_count++;
            if (left_trend_count >= PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT) {
                break;
            }
        } else {
            // Not positive => inconsistency
            inconsistency_count_left++;
            if (inconsistency_count_left > PEAK_VERIFICATION_ALLOWABLE_INCONSISTENCY_COUNT) {
                break;
            }
        }
    }

    // If we used nearly all available points to accumulate left trend, mark truncation
    if ((result.peak_index - left_trend_count) <= 0) {
        result.is_truncated_left = true;
    }

    /*************************************************************************
     * Right-Side Check: look for consistent NEGATIVE gradients.
     * We'll move right from `peak_index` up to `length - 1` or until we exceed 
     * allowable inconsistencies.
     *************************************************************************/
    for (uint16_t i = peak_index; i < (length - 1); ++i) {
        double grad = second_order_gradients[i + 1];
        if (grad < 0) {
            right_trend_count++;
            if (right_trend_count >= PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT) {
                break;
            }
        } else {
            // Not negative => inconsistency
            inconsistency_count_right++;
            if (inconsistency_count_right > PEAK_VERIFICATION_ALLOWABLE_INCONSISTENCY_COUNT) {
                break;
            }
        }
    }

    // If we used nearly all available points to accumulate right trend, mark truncation
    if ((result.peak_index + right_trend_count) >= (start_index + length - 1)) {
        result.is_truncated_right = true;
    }
    
     /*************************************************************************
     * Calculate Distance to End of Window
     *************************************************************************/
    uint16_t window_end_index = start_index + length - 1;
    result.distance_to_end = (result.peak_index <= window_end_index) ? 
                              (window_end_index - result.peak_index) : 0;

    printf("Distance from verified peak to end of window: %u points\n", result.distance_to_end);


    /*************************************************************************
     * Final Decision
     *************************************************************************/
    if ((left_trend_count  >= PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT) &&
        (right_trend_count >= PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT)) {
        result.peak_found = true;
        printf("Peak at index %u verified successfully.\n", result.peak_index);
    } else {
        printf("Peak at index %u failed verification.\n", result.peak_index);
    }

    return result;
}


/**
 * @brief Finds and verifies peaks in a data range using second-order gradients and a loop-based detection.
 *
 * This function:
 *   1) Adjusts the start index to avoid underflow.
 *   2) Computes second-order gradients over the adjusted interval via `trackGradients`.
 *   3) Loops through the gradients to detect sign changes (+ to -).
 *   4) Immediately verifies each detected candidate using `verify_quadratic_peak`.
 *   5) Stops when a verified peak is found or if the loop ends.
 *
 * @param data            Array of data points (MqsRawDataPoint_t) to analyze.
 * @param dataLength      Total number of data points in the `data` array.
 * @param increaseStart   Start index of the significant increase interval leading to the peak.
 * @param increaseEnd     End index of the significant increase interval where the peak is detected.
 * @param degree          Degree of the polynomial used for regression analysis.
 * @return QuadraticPeakAnalysisResult containing verification results.
 */
QuadraticPeakAnalysisResult verifyPeakValidity(
    const MqsRawDataPoint_t *data,
    uint16_t dataLength,
    uint16_t increaseStart,
    uint16_t increaseEnd,
    uint8_t degree
) {
    // Prepare the default result
    QuadraticPeakAnalysisResult result = {
        .peak_found        = false,
        .peak_index        = 0,
        .is_truncated_left = false,
        .is_truncated_right= false
    };

    // 1) Adjust start index to avoid underflow
    uint16_t adjustedStartIndex = (increaseStart >= START_INDEX_OFFSET)
                                ? (uint16_t)(increaseStart - START_INDEX_OFFSET)
                                : 0;

    // 2) Prepare to compute second-order gradients
    GradientCalculationResult secondOrderGradients;
    memset(&secondOrderGradients, 0, sizeof(GradientCalculationResult));

    // 3) Compute second-order gradients via trackGradients
    trackGradients(
        data,                      // Now passing MqsRawDataPoint_t*
        dataLength,
        adjustedStartIndex,
        degree,
        calculate_second_order_gradient,
        &secondOrderGradients
    );

    // 4) Loop to detect sign changes in secondOrderGradients
    //    (+ to - indicates a potential peak)
    for (uint16_t i = 1; i < secondOrderGradients.size; ++i) {
        double prevGradient = secondOrderGradients.gradients[i - 1];
        double currGradient = secondOrderGradients.gradients[i];

        if (prevGradient > 0.0 && currGradient < 0.0) {
            // Candidate peak => verify immediately using verify_quadratic_peak
            result = verify_quadratic_peak(
                secondOrderGradients.gradients,   // pass the array of double gradients
                secondOrderGradients.size,        // length of that array
                i,                                // relative peak index in the gradient array
                adjustedStartIndex                // offset for absolute indexing
            );

            if (result.peak_found) {
                DEBUG_PRINT_1("Verified peak found at absolute index %u.\n", result.peak_index);
                break;  // Stop once we have a verified peak
            } else {
                DEBUG_PRINT_1("Candidate at relative index %u (abs %u) failed verification.\n",
                       i, (uint16_t)(adjustedStartIndex + i));
            }
        }
    }

    // 5) If no verified peak found, result.peak_found remains false
    if (!result.peak_found) {
        printf("No verified peak found in the specified window.\n");
    }

    return result;
}


/**
 * @brief Calculates the required window shift to center the data around a peak based on trend durations.
 *
 * This function determines the direction and magnitude needed to adjust the data window
 * to better center it around a detected peak. It compares the lengths of consistent
 * increasing and decreasing trends and calculates the necessary shift based on their
 * difference. The window should be moved towards the side with the longer trend to achieve centering.
 *
 * **Logic:**
 * - If the increasing trend is longer than the decreasing trend, shift the window to the right.
 * - If the decreasing trend is longer than the increasing trend, shift the window to the left.
 * - If both trends are equal in length, the window is considered centered on the peak.
 *
 * **Rationale:**
 * - Shifting the window towards the side with the longer trend helps in aligning the window
 *   more symmetrically around the peak, ensuring that the analysis focuses on the most significant
 *   portion of the peak.
 *
 * @param trendResult Pointer to the `GradientTrendResultAbsolute` containing trend information.
 *
 * @return PeakCenteringResult  
 * A structure containing:
 * - `moveDirection`: Direction to move the window (`MOVE_LEFT`, `MOVE_RIGHT`, `ON_THE_PEAK`, `NO_TREND`).
 * - `moveAmount`: Amount by which to shift the window.
 * - `isCentered`: Indicates if the window is already centered on the peak.
 */
PeakCenteringResult calculatePeakCentering(const GradientTrendResultAbsolute *trendResult) {
    PeakCenteringResult centerResult;
    centerResult.isCentered = false;
    centerResult.moveAmount = 0;
    centerResult.moveDirection = NO_TREND;  // Default value

    // Constants

    // Calculate the lengths of the increasing and decreasing trends
    uint16_t increaseDuration = 0;
    uint16_t decreaseDuration = 0;

    if (trendResult->absoluteIncrease.startIndex != UINT16_MAX && 
        trendResult->absoluteIncrease.endIndex != UINT16_MAX) {
        increaseDuration = trendResult->absoluteIncrease.endIndex - trendResult->absoluteIncrease.startIndex + 1;
    }

    if (trendResult->absoluteDecrease.startIndex != UINT16_MAX && 
        trendResult->absoluteDecrease.endIndex != UINT16_MAX) {
        decreaseDuration = trendResult->absoluteDecrease.endIndex - trendResult->absoluteDecrease.startIndex + 1;
    }

    // Determine if the window is already centered
    if (increaseDuration == decreaseDuration && increaseDuration > 0) {
        centerResult.isCentered = true;
        centerResult.moveDirection = ON_THE_PEAK;
        centerResult.moveAmount = 0;
        DEBUG_PRINT_1("Window is already centered on the peak.\n");
        return centerResult;
    }

    // Determine shift amount and direction based on trend durations
    int shiftAmount = 0;
    if (increaseDuration > decreaseDuration) {
        shiftAmount = (increaseDuration - decreaseDuration) / CENTERING_RATIO;
        centerResult.moveDirection = MOVE_RIGHT;
        centerResult.moveAmount = shiftAmount;
        DEBUG_PRINT_1("Increase duration (%u) > decrease duration (%u). Moving right by %d.\n",
                      increaseDuration, decreaseDuration, shiftAmount);
    }
    else if (decreaseDuration > increaseDuration) {
        shiftAmount = (decreaseDuration - increaseDuration) / CENTERING_RATIO;
        centerResult.moveDirection = MOVE_LEFT;
        centerResult.moveAmount = shiftAmount;
        DEBUG_PRINT_1("Decrease duration (%u) > increase duration (%u). Moving left by %d.\n",
                      decreaseDuration, increaseDuration, shiftAmount);
    }
    else {
        // If both durations are zero or unequal but zero, consider it centered
        centerResult.isCentered = true;
        centerResult.moveDirection = ON_THE_PEAK;
        centerResult.moveAmount = 0;
        DEBUG_PRINT_1("Trend durations are equal or zero. Window is centered.\n");
    }

    return centerResult;
}

/**
 * @brief Identifies gradient trends and calculates peak centering parameters.
 *
 * This function performs the following steps:
 * 1. Identifies first-order gradient trends within the specified analysis window.
 * 2. Calculates the required window shift to center around the detected peak based on trend durations.
 *
 * **Note:** This function does not modify the `startIndex`. The caller is responsible for adjusting the window position.
 *
 * @param data          Array of MqsRawDataPoint_t (phase angle data) to analyze.
 * @param startIndex    Starting index of the analysis window.
 * @param analysisLength Length of the analysis window.
 * @param degree        Degree of the polynomial used for regression.
 *
 * @return PeakCenteringResult  
 * A structure containing:
 *  - `moveDirection`: Direction to move the window (`MOVE_LEFT`, `MOVE_RIGHT`, `ON_THE_PEAK`, `NO_TREND`).
 *  - `moveAmount`: Amount by which to shift the window.
 *  - `isCentered`: Indicates if the window is already centered on the peak.
 */
PeakCenteringResult identifyAndCalculateCentering(
    const MqsRawDataPoint_t *data,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree
) {
    // 1. Identify first-order gradient trends
    GradientTrendResultAbsolute trendResultFirstOrder;
    GradientCalculationResult   gradCalcResultFirstOrderInc;
    GradientCalculationResult   gradCalcResultFirstOrderDec;

    // Initialize structures
    initializeGradientTrendResultAbsolute(&trendResultFirstOrder);
    memset(&gradCalcResultFirstOrderInc, 0, sizeof(GradientCalculationResult));
    memset(&gradCalcResultFirstOrderDec, 0, sizeof(GradientCalculationResult));

    // Use identifyGradientTrends on the MqsRawDataPoint_t array
    identifyGradientTrends(
        data,
        startIndex,
        analysisLength,
        degree,
        GRADIENT_ORDER_FIRST,
        TREND_TYPE_GRADIENT_TRENDS,
        &trendResultFirstOrder,
        &gradCalcResultFirstOrderInc,
        &gradCalcResultFirstOrderDec
    );

    // 2. Calculate peak centering
    PeakCenteringResult centerResult = calculatePeakCentering(&trendResultFirstOrder);

    // Debug info
    if (centerResult.isCentered) {
        DEBUG_PRINT_1("Window is already centered on the peak at startIndex=%u.\n", startIndex);
    } else {
        // Output the move direction and amount
        printf("Move Direction: ");
        switch (centerResult.moveDirection) {
            case MOVE_RIGHT:  printf("MOVE_RIGHT\n");  break;
            case MOVE_LEFT:   printf("MOVE_LEFT\n");   break;
            case ON_THE_PEAK: printf("ON_THE_PEAK\n"); break;
            case NO_TREND:    printf("NO_TREND\n");    break;
            default:          printf("UNKNOWN\n");     break;
        }
        printf("Move Amount: %u\n", centerResult.moveAmount);
    }

    return centerResult;
}

/**
 * @brief Checks for a continuous negative trend in the second-order gradients,
 *        ensuring the trend starts after a specified minimum ratio of the window
 *        and reaches a defined negative gradient threshold.
 *
 * This function computes second-order gradients for a dataset and determines whether
 * there is a consecutive negative trend of at least 5 values that starts after a
 * specified portion (`min_ratio`) of the window and reaches a gradient of `-1.0` or below.
 * The results are stored in the provided `GradientCalculationResult` structure.
 *
 * @param data          Array of MqsRawDataPoint_t to analyze (reading .phaseAngle from each).
 * @param length        The number of data points in the `data` array.
 * @param start_index   The starting index in the `data` array to begin analysis.
 * @param min_ratio     The minimum ratio (e.g., 0.3 for 30%) indicating where the negative
 *                      trend should start within the window.
 * @param result        Pointer to a `GradientCalculationResult` structure to store the results.
 * @return `true` if a valid negative trend is found, starts after `min_ratio` of the window,
 *         and reaches the threshold of `-1.0` or below; `false` otherwise.
 */
bool peakAnalysis_checkPeakCentering(
    const MqsRawDataPoint_t *data,
    uint16_t length,
    uint16_t start_index,
    double min_ratio,
    GradientCalculationResult *result
) {
    DEBUG_PRINT_3("Entering peakAnalysis_checkPeakCentering with start_index=%u, length=%u, min_ratio=%.2f\n",
                  start_index, length, min_ratio);

    // Initialize the result structure
    result->size = 0;
    result->startIndex = 0;
    result->endIndex = 0;
    result->median = 0.0;  
    result->mad = 0.0;     

    // 1) Track gradients using second-order polynomial regression
    //    but we'll actually calculate first-order? 
    //    (As in your original code, it's "degree=3" plus "calculate_first_order_gradient".)
    //    Adjust if you truly want second-order logic.
    trackGradients(
        data,                        // MqsRawDataPoint_t array
        length,
        start_index,
        3,  // Degree = 3 for second-order gradients
        calculate_first_order_gradient, // Possibly rename if you truly want second-order gradient
        result
    );

    DEBUG_PRINT_3("Collected %u gradients.\n", result->size);

    // 2) Ensure we have enough gradients to detect a negative trend
    if (result->size < PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT) { 
        DEBUG_PRINT_1("Not enough gradients to determine a negative trend.\n");
        return false;
    }

    // 3) Calculate the minimum acceptable index based on min_ratio
    double min_acceptable_index_f = min_ratio * (double)result->size;
    uint16_t min_acceptable_index = (uint16_t)ceil(min_acceptable_index_f); 
    DEBUG_PRINT_3("Minimum acceptable trend start index: %u (min_ratio=%.2f)\n", 
                  min_acceptable_index, min_ratio);

    // 4) Search for a continuous negative trend of >= 5 consecutive negative gradients
    //    that begins after min_acceptable_index and also hits <= -1.0
    uint16_t consecutive_negatives = 0;
    uint16_t trend_start_index     = 0;

    for (uint16_t i = 0; i < result->size; ++i) {
        double gradValue = result->gradients[i];
        DEBUG_PRINT_2("Inspecting gradient[%u] = %.2f\n", i, gradValue);

        if (gradValue < 0.0) {
            // Potential part of a negative trend
            if (consecutive_negatives == 0) {
                // Mark new potential trend start
                trend_start_index = i;
                DEBUG_PRINT_2("Potential new negative trend starting at index %u.\n", trend_start_index);
            }
            consecutive_negatives++;
            DEBUG_PRINT_2("Consecutive negatives count: %u\n", consecutive_negatives);

            // Check if we have at least 5 consecutive negative gradients
            if (consecutive_negatives >= PEAK_VERIFICATION_MIN_CONSISTENT_TREND_COUNT) {
                DEBUG_PRINT_2("Found %u consecutive negative gradients starting at index %u.\n",
                              consecutive_negatives, trend_start_index);

                // Verify if the trend starts after the min_ratio threshold
                if (trend_start_index >= min_acceptable_index) {
                    DEBUG_PRINT_2("Trend starts after min_ratio threshold. Verifying further.\n");

                    // Check if any gradient in [trend_start_index..i] <= -1.0
                    bool threshold_reached = false;
				
#if DEBUG_LEVEL > 0
										uint16_t threshold_index = 0; // Declare only if DEBUG_LEVEL > 0
#endif

										for (uint16_t j = trend_start_index; j <= i; ++j) {
												double gVal = result->gradients[j];
												DEBUG_PRINT_3("Checking gradient[%u] = %.2f against threshold %.2f\n", 
																			j, gVal, NEGATIVE_GRADIENT_THRESHOLD);
												if (gVal <= NEGATIVE_GRADIENT_THRESHOLD) {
														threshold_reached = true;
#if DEBUG_LEVEL > 0
														threshold_index = j; // Initialize only if DEBUG_LEVEL > 0
														DEBUG_PRINT_2("Threshold %.2f reached at index %u.\n", NEGATIVE_GRADIENT_THRESHOLD, threshold_index);
#endif
														break;
												}
										}
                    if (threshold_reached) {
                        // Mark the start/end indices in the original data
                        result->startIndex = (uint16_t)(start_index + trend_start_index);
                        result->endIndex   = (uint16_t)(start_index + i);
                        printf("Found valid negative trend: startIndex=%u, endIndex=%u, gradient=%.2f.\n", 
                               result->startIndex, result->endIndex, gradValue);
                        return true;
                    } else {
                        DEBUG_PRINT_2("Negative trend starts after min_ratio but does not reach threshold of %.2f.\n",
                                      NEGATIVE_GRADIENT_THRESHOLD);
                        // Reset and continue searching
                        consecutive_negatives = 0;
                    }
                } else {
                    DEBUG_PRINT_2("Negative trend starts at index %u, which is before min_acceptable_index %u.\n",
                                  trend_start_index, min_acceptable_index);
                    consecutive_negatives = 0;
                }
            }
        } else {
            // If the current gradient is >= 0, the negative sequence is broken
            if (consecutive_negatives > 0) {
                DEBUG_PRINT_3("Negative trend broken at index %u after %u consecutive negatives.\n",
                              i, consecutive_negatives);
            }
            consecutive_negatives = 0;
        }
    }

    // If we finish the loop, we didn't find a valid negative trend
    DEBUG_PRINT_1("No valid negative trend found that starts after the minimum ratio and reaches threshold.\n");
    return false;
}

// 4 functions exist for raw data point conversion. 

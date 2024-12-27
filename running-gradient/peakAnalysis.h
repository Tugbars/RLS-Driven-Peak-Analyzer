/**
 * @file peakSignificanceAnalysis.h
 * @brief Header file for peak significance analysis functions.
 *
 * This header defines the structures and function prototypes for analyzing the significance
 * of detected peaks within gradient trends.
 */

#ifndef PEAK_SIGNIFICANCE_ANALYSIS_H
#define PEAK_SIGNIFICANCE_ANALYSIS_H

#include "trend_detection.h"  // To access GradientTrendResult and related structures
#include <stdint.h>
#include <stdbool.h>

#include "mqs_def.h"

/**
 * @brief Structure to store the results of peak centering calculations.
 */
typedef struct {
    MoveDirection moveDirection;  /**< Direction to move the window */
    uint16_t moveAmount;          /**< Amount by which to shift the window */
    bool isCentered;              /**< Indicates if the window is already centered on the peak */
} PeakCenteringResult;

/**
 * @brief Structure to hold peak analysis results.
 */
typedef struct {
    bool isOnPeak;                 /**< Indicates if the current window is on a peak */
    bool isSignificantPeak;        /**< Indicates if a significant peak was detected */
    bool isWindowCentered;         /**< Indicates if the window is centered enough based on second-order gradients */
    bool isValidPeak;              /**< Indicates if the peak is valid based on second-order gradient check */
    PeakPosition moveDirection;   /**< Direction to move the window if centering is needed */
} PeakAnalysisResult;


/**
 * @brief Structure to hold the result of quadratic peak verification.
 */
typedef struct {
    bool peak_found;            /**< Indicates if a peak was found and verified */
    uint16_t peak_index;        /**< Index of the detected peak in the dataset */
    bool is_truncated_left;     /**< Indicates if the left side of the peak is truncated */
    bool is_truncated_right;    /**< Indicates if the right side of the peak is truncated */
} QuadraticPeakAnalysisResult;

/**
 * @brief Performs peak analysis by identifying gradient trends and analyzing them.
 *
 * @param values             Array of data points to analyze.
 * @param dataLength         The total length of the data array.
 * @param startIndex         The starting index in the values array.
 * @param analysisLength     The number of data points to include in the analysis window.
 * @param degree             The degree of the polynomial used for regression.
 * @param gradientOrder      The order of the gradient (first-order or second-order).
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
);


/**
 * @brief Validates a detected peak by analyzing its curvature using second-order gradients.
 *
 * The `verifyPeakValidity` function confirms the authenticity of a detected peak within a dataset
 * by examining the second-order gradients (curvature) around the peak. This verification ensures
 * that the peak exhibits the necessary concave-down characteristics on both sides, distinguishing
 * genuine peaks from random fluctuations or noise.

 * **Rationale:**
 * - **Second-Order Gradients:** Provide insight into the curvature of the data, ensuring the peak has a
 *   concave-down shape characteristic of true peaks.
 * - **Trend Consistency:** Verifies that the peak is part of a broader increasing and then decreasing
 *   trend, reducing the likelihood of false positives caused by isolated data points.
 *
 * @param values           Array of data points to analyze.
 * @param dataLength       Total number of data points in the `values` array.
 * @param increaseStart    Start index of the significant increase interval leading to the peak.
 * @param increaseEnd      End index of the significant increase interval where the peak is detected.
 * @param degree           Degree of the polynomial used for regression analysis.
 */
QuadraticPeakAnalysisResult verifyPeakValidity(
    const MqsRawDataPoint_t *data,
    uint16_t dataLength,
    uint16_t increaseStart,
    uint16_t increaseEnd,
    uint8_t degree
);

/**
 * @brief Calculates the required window shift to center the data around a peak based on trend durations.
 *
 * @param trendResult Pointer to the `GradientTrendResultAbsolute` containing trend information.
 *
 * @return PeakCenteringResult  
 * A structure containing:
 * - `moveDirection`: Direction to move the window (`MOVE_LEFT`, `MOVE_RIGHT`, `ON_THE_PEAK`, `NO_TREND`).
 * - `moveAmount`: Amount by which to shift the window.
 * - `isCentered`: Indicates if the window is already centered on the peak.
 */
PeakCenteringResult calculatePeakCentering(const GradientTrendResultAbsolute *trendResult);

/**
 * @brief Identifies gradient trends and calculates peak centering parameters.
 *
 * This function performs the following steps:
 * 1. Identifies first-order gradient trends within the specified analysis window.
 * 2. Calculates the required window shift to center around the detected peak based on trend durations.
 * 
 * @param phaseAngles    Array of phase angle data points.
 * @param startIndex     Starting index of the analysis window.
 * @param analysisLength Length of the analysis window.
 * @param degree         Degree of the polynomial used for regression.
 */
PeakCenteringResult identifyAndCalculateCentering(
    const MqsRawDataPoint_t *data,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree
);

/**
 * @brief Centers the data window around a detected peak by calculating shift parameters.
 *
 * This function performs the following:
 * 1. Identifies gradient trends within the current analysis window.
 * 2. Calculates the direction and amount to shift the window to center around the peak.
 *
 * @param phaseAngles    Array of phase angle data points.
 * @param dataLength     Total number of data points in `phaseAngles`.
 * @param startIndex     Pointer to the starting index of the analysis window.
 * @param analysisLength Length of the analysis window.
 * @param degree         Degree of the polynomial used for regression.
 */
PeakCenteringResult centerWindowOnPeak(
    double *phaseAngles,
    uint16_t dataLength,
    uint16_t *startIndex,
    uint16_t analysisLength,
    uint8_t degree
);

/**
 * @brief Checks for a continuous negative trend in the second-order gradients,
 *        ensuring the trend starts after a specified minimum ratio of the window.
 *
 * This function computes second-order gradients for a dataset and determines whether
 * there is a consecutive negative trend of at least 5 values that starts after a
 * specified portion (`min_ratio`) of the window. The results are stored in the provided
 * `GradientCalculationResult` structure.
 *
 * @param measurements  Array of measurement values to analyze.
 * @param length        The number of data points in the `measurements` array.
 * @param start_index   The starting index in the `measurements` array to begin analysis.
 * @param min_ratio     The minimum ratio (e.g., 0.3 for 30%) indicating where the negative
 *                      trend should start within the window.
 * @param result        Pointer to a `GradientCalculationResult` structure to store the results.
 * @return `true` if a valid negative trend is found and starts after `min_ratio` of the window,
 *         `false` otherwise.
 */
bool peakAnalysis_checkPeakCentering(
    const MqsRawDataPoint_t *data,
    uint16_t length,
    uint16_t start_index,
    double min_ratio,
    GradientCalculationResult *result
);

#endif // PEAK_SIGNIFICANCE_ANALYSIS_H

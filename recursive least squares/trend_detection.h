/**
 * @file trend_detection.h
 * @brief Defines trend detection structures and functions for identifying consistent increasing and decreasing trends
 *        using gradients calculated via Recursive Least Squares (RLS) polynomial regression.
 *
 * This header provides definitions for structures and function prototypes used in trend detection.
 * It facilitates the identification and analysis of consistent trends within a dataset of gradient values.
 * 
 * The trend detection mechanisms are essential for applications like signal processing, financial data analysis,
 * or any time-series data analysis where discerning underlying trends is crucial.
 *
 * The functions utilize Recursive Least Squares (RLS) regression from `rls_polynomial_regression.h`
 * to calculate first-order and second-order gradients over a sliding window, offering a robust
 * method for online trend analysis.
 */

#ifndef TREND_DETECTION_H
#define TREND_DETECTION_H

#include <stdint.h>
#include <stdbool.h>
#include "rls_polynomial_regression.h" // Ensure this header file exists and is correctly implemented
#include "mqs_def.h"

/** 
 * @def MAX_WINDOW_SIZE
 * @brief Maximum size of the window for trend detection.
 */
#define MAX_WINDOW_SIZE 30 // Adjust as per application requirements

/** 
 * @def MAX_TREND_SIZE
 * @brief Maximum size of a detected trend.
 */
#define MAX_TREND_SIZE 30 // Adjust as necessary

/**
 * @brief Enumeration for different types of trend detection.
 */
typedef enum {
    TREND_TYPE_ABSOLUTE,           /**< Detects absolute increasing and decreasing trends based on median. */
    TREND_TYPE_SIGNIFICANT,        /**< Detects significant increasing and decreasing trends based on median and MAD. */
    TREND_TYPE_GRADIENT_TRENDS     /**< Custom trend detection type for gradient-based trends. */
} TrendDetectionType;

typedef enum {
    GRADIENT_ORDER_FIRST,
    GRADIENT_ORDER_SECOND
} GradientOrder;


/**
 * @enum PeakPosition
 * @brief An enumeration to represent the position relative to a peak.
 *
 * This enumeration is used to indicate whether a detected trend is on the left side, right side, or undecided relative to a peak.
 */
typedef enum {
    LEFT_SIDE = -1,                       /**< Indicates that the trend is on the left side of a peak. */
    RIGHT_SIDE = 1,                       /**< Indicates that the trend is on the right side of a peak. */
    UNDECIDED = 0,                        /**< Indicates that the trend's position relative to a peak is undecided. */
    ON_PEAK = 2,                          /**< Indicates that the trend is on the peak. */
    NEGATIVE_UNDECIDED = 3                /**< Indicates that both sides have negative trends with no significant difference. */
} PeakPosition;

/**
 * @brief Enumeration for move directions based on trend analysis.
 */
typedef enum {
    MOVE_RIGHT,    /**< Indicates shifting the window to the right towards increasing trends. */
    MOVE_LEFT,     /**< Indicates shifting the window to the left towards decreasing trends. */
    ON_THE_PEAK,   /**< Indicates that the current window is at a significant peak. */
    NO_TREND       /**< Indicates no dominant trend detected. */
} MoveDirection;

/**
 * @brief Structure to hold the indices and values for a gradient trend.
 */
typedef struct {
    uint16_t startIndex;           /**< Start index of the trend */
    uint16_t endIndex;             /**< End index of the trend */
    double maxValue;               /**< Maximum (for increase) or minimum (for decrease) gradient value within the trend */
    uint16_t maxValueIndex;        /**< Index of the maximum or minimum gradient value within the trend */
    double sumGradients;           /**< Sum of gradients within the trend */
} GradientTrendIndices;

/**
 * @brief Base structure for gradient trend results.
 *
 * This structure provides an abstraction over different trend result structures by holding
 * pointers to their respective fields. It allows functions to operate on trend results
 * generically without knowing the specific type of trend result structure.
 */
typedef struct {
    GradientTrendIndices *increase;      /**< Pointer to the increase trend information */
    GradientTrendIndices *decrease;      /**< Pointer to the decrease trend information */
    PeakPosition *moveDirection;         /**< Pointer to the move direction */
    bool *isSignificantPeak;             /**< Pointer to the significant peak flag */
    int16_t *moveAmount;                 /**< Pointer to the move amount */
} GradientTrendResultBase;

/**
 * @brief Structure to hold the result of absolute gradient trend detection.
 */
typedef struct {
    GradientTrendIndices absoluteIncrease;    /**< Information about the longest consistent increasing trend */
    GradientTrendIndices absoluteDecrease;    /**< Information about the longest consistent decreasing trend */
    PeakPosition moveDirection;               /**< Move direction based on the absolute trend */
    bool isSignificantPeak;                   /**< Indicates if a significant peak was detected */
    int16_t moveAmountAbsolute;               /**< Amount to move based on absolute trend */
    GradientTrendResultBase base;             /**< Base structure pointing to the above fields */
} GradientTrendResultAbsolute;

/**
 * @brief Structure to hold the result of significant gradient trend detection.
 */
typedef struct {
    GradientTrendIndices significantIncrease; /**< Information about the longest consistent significant increasing trend */
    GradientTrendIndices significantDecrease; /**< Information about the longest consistent significant decreasing trend */
    PeakPosition moveDirection;               /**< Move direction based on the significant trend */
    bool isSignificantPeak;                   /**< Indicates if a significant peak was detected */
    int16_t moveAmountSignificant;            /**< Amount to move based on significant trend */
    GradientTrendResultBase base;             /**< Base structure pointing to the above fields */
} GradientTrendResultSignificant;

/**
 * @brief Initializes the GradientTrendResultAbsolute structure.
 *
 * This function resets the fields of the GradientTrendResultAbsolute structure to their default values.
 * It also sets up the base pointers to the respective fields, allowing for generic access in functions
 * that operate on the base structure.
 *
 * @param result Pointer to the GradientTrendResultAbsolute structure to initialize.
 */
void initializeGradientTrendResultAbsolute(GradientTrendResultAbsolute *result);

/**
 * @brief Initializes the GradientTrendResultSignificant structure.
 *
 * Similar to `initializeGradientTrendResultAbsolute`, this function resets the fields of the
 * GradientTrendResultSignificant structure and sets up the base pointers.
 *
 * @param result Pointer to the GradientTrendResultSignificant structure to initialize.
 */
void initializeGradientTrendResultSignificant(GradientTrendResultSignificant *result);

/**
 * @brief Detects trends (absolute or significant) and calculates move amounts using separate statistics for increases and decreases.
 *
 * This function analyzes the gradient calculation results to detect increasing and decreasing trends.
 * It uses separate medians and MADs for increasing and decreasing trends to determine significant trends.
 * The function populates a trend result structure with the detected trends and move amounts.
 *
 * @param gradCalcResultInc Pointer to the gradient calculation results for increases.
 * @param gradCalcResultDec Pointer to the gradient calculation results for decreases.
 * @param trendResultBase Pointer to the base of the trend result structure.
 * @param trendType Specifies whether to detect absolute or significant trends.
 */
void detectTrends(
    const GradientCalculationResult *gradCalcResultInc,
    const GradientCalculationResult *gradCalcResultDec,
    GradientTrendResultBase *trendResultBase,
    TrendDetectionType trendType
);


/**
 * @brief Identifies gradient trends within a window using RLS regression and MAD-based statistics.
 *
 * This function calculates gradients using the provided `calculate_gradient` function, identifies the
 * longest consistent increasing and decreasing trends using `findConsistentIncreaseInternal` and
 * `findConsistentDecreaseInternal`, computes separate statistical measures (median and MAD)
 * for each trend, and stores the results in the provided `gradCalcResultInc` and `gradCalcResultDec` structures.
 *
 * @param values            Array of data points to analyze.
 * @param startIndex        The starting index in the values array.
 * @param analysisLength    The number of data points to include in the gradient analysis window.
 * @param degree            The degree of the polynomial used for regression.
 * @param calculate_gradient Function pointer to the gradient calculation function.
 * @param gradCalcResultInc Pointer to store the gradient calculation results for increases.
 * @param gradCalcResultDec Pointer to store the gradient calculation results for decreases.
 *
 * @return void
 */
PeakPosition identifyTrends(
    const MqsRawDataPoint_t *data,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree,
    double (*calculate_gradient)(const RegressionState *, double),
    GradientCalculationResult *gradCalcResultInc,
    GradientCalculationResult *gradCalcResultDec
);

/**
 * @brief Identifies gradient trends within a specified window.
 *
 * This function detects consistent increasing and decreasing trends based on first-order or second-order gradients
 * calculated using RLS regression. It delegates the detection to specific functions based on the
 * trend detection type (absolute, significant, or gradient trends).
 *
 * @param values             Array of data points to analyze.
 * @param startIndex         The starting index in the values array.
 * @param analysisLength     The number of data points to include in the gradient analysis window.
 * @param degree             The degree of the polynomial used for regression.
 * @param gradientOrder      The order of the gradient (first-order or second-order).
 * @param trendType          The type of trend detection to perform.
 * @param trendResult        Pointer to store the trend detection result.
 * @param gradCalcResultInc  Pointer to store the gradient calculation results for increases.
 * @param gradCalcResultDec  Pointer to store the gradient calculation results for decreases.
 */
void identifyGradientTrends(
    const MqsRawDataPoint_t *data,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree,
    GradientOrder gradientOrder,
    TrendDetectionType trendType,
    void *trendResult,
    GradientCalculationResult *gradCalcResultInc,
    GradientCalculationResult *gradCalcResultDec
);

/**
 * @brief Wrapper function that identifies trends and detects the move direction.
 *
 * This function calculates the gradients using RLS polynomial regression by calling `identifyTrends`,
 * and then it uses those gradients to determine the move direction by calling `detectGradientTrends`.
 * It encapsulates the entire process of gradient calculation and trend detection, simplifying the usage.
 *
 * @param values             Array of data points to analyze.
 * @param startIndex         The starting index in the values array.
 * @param analysisLength     The number of data points to include in the gradient analysis window.
 * @param degree             The degree of the polynomial used for regression.
 * @param trendType          The type of trend detection to perform.
 * @param trendResult        Pointer to store the trend detection result.
 * @param calculate_gradient Function pointer to the gradient calculation function.
 *
 * @return MoveDirection The determined move direction based on trend analysis.
 */
MoveDirection identifyAndDetectTrends(
    const double *values,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree,
    GradientOrder gradientOrder,
    TrendDetectionType trendType,
    void *trendResult
);

/**
 * @brief Identifies the move direction based on gradient trends by splitting the gradient array into halves.
 *
 * This function calculates the gradients using `trackGradients`, splits the gradient array into two halves,
 * assesses each part to determine the direction towards the impedance peak, and returns the corresponding move direction.
 *
 * @param values          Array of data points to analyze.
 * @param startIndex      The starting index in the values array.
 * @param analysisLength  The number of data points to include in the gradient analysis window.
 * @param degree          The degree of the polynomial used for regression.
 * @param gradientOrder   The order of the gradient (first-order or second-order).
 * @param trendType       The type of trend detection to perform (currently unused in this implementation).
 * @param trendResult     Pointer to store the trend detection result (optional; not used here).
 *
 * @return MoveDirection The determined move direction based on trend analysis.
 */
MoveDirection identifySplitMoveDirection(
    const MqsRawDataPoint_t *data,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree,
    GradientOrder gradientOrder
);

/**
 * @brief Adjusts the window start index based on the move amount.
 *
 * This function calculates a new start index for the analysis window by adding the move amount
 * to the current start index. It ensures that the new start index stays within the bounds
 * of the dataset, preventing out-of-bounds access.
 *
 * @param currentStartIndex The current start index of the window.
 * @param moveAmount The amount to move the window (positive or negative).
 * @param dataLength The total length of the data array.
 * @param windowSize The size of the window.
 * @return uint16_t The new start index after adjustment.
 */
uint16_t adjustWindowStartIndex(
    uint16_t currentStartIndex,
    int16_t moveAmount,
    uint16_t dataLength,
    uint16_t windowSize
);

void printTrendDetails(const char* trendName, const GradientTrendIndices* trendIndices);

#endif // TREND_DETECTION_H



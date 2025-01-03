/**
 * @file trend_detection.c
 * @brief Implements trend detection functions for identifying consistent increasing and decreasing trends
 *        using gradients calculated via Recursive Least Squares (RLS) polynomial regression.
 *
 * This file provides functions to identify and analyze consistent trends within a dataset of gradient values.
 * It separately tracks the maximum and minimum gradients within increasing and decreasing trends,
 * respectively, and accounts for minor fluctuations using configurable parameters. 
 * 
 * The trend detection mechanisms are essential for applications like signal processing, financial data analysis,
 * or any time-series data analysis where discerning underlying trends is crucial.
 *
 * The functions utilize Recursive Least Squares (RLS) regression from `rls_polynomial_regression.h`
 * to calculate first-order and second-order gradients over a sliding window, offering a robust
 * method for online trend analysis.Calculating
 */

#include "trend_detection.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>  // For va_list, va_start, va_end
#include "statistics.h"

// Preprocessor constants for trend detection configuration

/** 
 * @def SIGNIFICANCE_FACTOR
 * @brief Factor multiplied by MAD to determine significance threshold.
 *
 * This factor adjusts the sensitivity of trend detection. A higher value ensures that only larger deviations
 * from the median are considered significant.
 */
#define SIGNIFICANCE_FACTOR 0.3  // Adjust this factor as needed

/** 
 * @def MAX_ALLOWED_DECREASES
 * @brief Maximum number of consecutive minor decreases allowed within an increasing trend before considering the trend broken.
 *
 * This parameter allows the trend detection algorithm to tolerate small fluctuations within an increasing trend.
 * If the number of consecutive decreases exceeds this value, the trend is considered to have ended.
 */
#define MAX_ALLOWED_DECREASES 1

/** 
 * @def MAX_ALLOWED_INCREASES
 * @brief Maximum number of consecutive minor increases allowed within a decreasing trend before considering the trend broken.
 *
 * Similar to `MAX_ALLOWED_DECREASES`, this parameter allows small increases within a decreasing trend.
 */
#define MAX_ALLOWED_INCREASES 1

/** 
 * @def TREND_THRESHOLD
 * @brief Threshold for determining significant trend counts.
 */
#define TREND_THRESHOLD (RLS_WINDOW / 4) // For example, 30 / 4 = 7

/** 
 * @def MIN_GRADIENT_THRESHOLD
 * @brief Minimum gradient threshold for undecided move direction.
 */
#define MIN_GRADIENT_THRESHOLD -0.25

/** 
 * @def MAX_GRADIENT_THRESHOLD
 * @brief Maximum gradient threshold for undecided move direction.
 */
#define MAX_GRADIENT_THRESHOLD 0.25

/** 
 * @def MIN_TREND_COUNT_FOR_PEAK
 * @brief Minimum trend count required to consider being on the peak.
 */
#define MIN_TREND_COUNT_FOR_PEAK 3

// Define debug levels
#define DEBUG_LEVEL 0  // Set to 0, 1, 2, or 3 to enable different levels of debugging

// Define separate debug print macros based on DEBUG_LEVEL

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

// Function to print trend details
void printTrendDetails(const char* trendName, const GradientTrendIndices* trendIndices) {
    if (trendIndices->startIndex != UINT16_MAX && trendIndices->endIndex != UINT16_MAX) {
        printf("%s Trend: ", trendName);
        printf("Start Index: %u ", trendIndices->startIndex);
        printf("End Index: %u ", trendIndices->endIndex);
        printf("Max Value: %.6f at Index %u\n", trendIndices->maxValue, trendIndices->maxValueIndex);
        //printf("Sum of Gradients: %.6f\n", trendIndices->sumGradients);
    } else {
        //printf("No %s Trend Detected.\n", trendName);
    }
}

/**
 * @brief Initializes the GradientTrendResultAbsolute structure.
 *
 * This function resets the fields of the GradientTrendResultAbsolute structure to their default values.
 * It also sets up the base pointers to the respective fields, allowing for generic access in functions
 * that operate on the base structure.
 *
 * @param result Pointer to the GradientTrendResultAbsolute structure to initialize.
 */
void initializeGradientTrendResultAbsolute(GradientTrendResultAbsolute *result) {
    // Initialize fields
    memset(result, 0, sizeof(GradientTrendResultAbsolute));
    result->absoluteIncrease.startIndex = UINT16_MAX;
    result->absoluteIncrease.endIndex = UINT16_MAX;
    result->absoluteIncrease.maxValueIndex = UINT16_MAX;

    result->absoluteDecrease.startIndex = UINT16_MAX;
    result->absoluteDecrease.endIndex = UINT16_MAX;
    result->absoluteDecrease.maxValueIndex = UINT16_MAX;

    result->moveDirection = UNDECIDED;
    result->isSignificantPeak = false;
    result->moveAmountAbsolute = 0;

    // Set up base pointers
    result->base.increase = &result->absoluteIncrease;
    result->base.decrease = &result->absoluteDecrease;
    result->base.moveDirection = &result->moveDirection;
    result->base.isSignificantPeak = &result->isSignificantPeak;
    result->base.moveAmount = &result->moveAmountAbsolute;
}


/**
 * @brief Initializes the GradientTrendResultSignificant structure.
 *
 * Similar to `initializeGradientTrendResultAbsolute`, this function resets the fields of the
 * GradientTrendResultSignificant structure and sets up the base pointers.
 *
 * @param result Pointer to the GradientTrendResultSignificant structure to initialize.
 */
void initializeGradientTrendResultSignificant(GradientTrendResultSignificant *result) {
    // Initialize fields
    memset(result, 0, sizeof(GradientTrendResultSignificant));
    result->significantIncrease.startIndex = UINT16_MAX;
    result->significantIncrease.endIndex = UINT16_MAX;
    result->significantIncrease.maxValueIndex = UINT16_MAX;

    result->significantDecrease.startIndex = UINT16_MAX;
    result->significantDecrease.endIndex = UINT16_MAX;
    result->significantDecrease.maxValueIndex = UINT16_MAX;

    result->moveDirection = UNDECIDED;
    result->isSignificantPeak = false;
    result->moveAmountSignificant = 0;

    // Set up base pointers
    result->base.increase = &result->significantIncrease;
    result->base.decrease = &result->significantDecrease;
    result->base.moveDirection = &result->moveDirection;
    result->base.isSignificantPeak = &result->isSignificantPeak;
    result->base.moveAmount = &result->moveAmountSignificant;
}


/**
 * @brief Finds the longest consistent increasing trend in the gradient array based on a specified threshold.
 *
 * This function scans through an array of gradient values to identify the longest consistent increasing trend
 * that exceeds a specified threshold. It allows a specified number of minor consecutive decreases to account
 * for small fluctuations. The function returns a `GradientTrendIndices` struct containing the start and end
 * indices of the longest consistent increase, the maximum gradient value within that trend,
 * and the sum of gradients within the trend.
 *
 * @param gradients   The array of gradient values to analyze.
 * @param windowSize  The number of elements in the gradient array to analyze.
 * @param threshold   The threshold above which gradients are considered significant.
 * @return GradientTrendIndices A structure containing the trend information.
 */
/**
 * @brief Finds the longest consistent increasing trend in the gradient array based on a specified threshold.
 *
 * This function scans through an array of gradient values to identify the longest consistent increasing trend
 * that exceeds a specified threshold. It allows a specified number of minor consecutive decreases to account
 * for small fluctuations. The function returns a `GradientTrendIndices` struct containing the start and end
 * indices of the longest consistent increase, the maximum gradient value within that trend,
 * and the sum of gradients within the trend.
 *
 * @param gradients   The array of gradient values to analyze.
 * @param windowSize  The number of elements in the gradient array to analyze.
 * @param threshold   The threshold above which gradients are considered significant.
 * @return GradientTrendIndices A structure containing the trend information.
 */
static GradientTrendIndices findConsistentIncreaseInternal(
    const double *gradients,
    uint16_t windowSize,
    double threshold
) {
    DEBUG_PRINT_3("Entering findConsistentIncreaseInternal\n");

    // **Initialize the structure to store the longest increasing trend found**
    GradientTrendIndices longestIncreaseInfo;
    memset(&longestIncreaseInfo, 0, sizeof(GradientTrendIndices));
    longestIncreaseInfo.startIndex = UINT16_MAX;      // Initialize to UINT16_MAX to indicate invalid index
    longestIncreaseInfo.endIndex = UINT16_MAX;
    longestIncreaseInfo.maxValueIndex = UINT16_MAX;

    // **Initialize the structure to store the current increasing trend being tracked**
    GradientTrendIndices currentIncreaseInfo;
    memset(&currentIncreaseInfo, 0, sizeof(GradientTrendIndices));
    currentIncreaseInfo.startIndex = UINT16_MAX;
    currentIncreaseInfo.endIndex = UINT16_MAX;
    currentIncreaseInfo.maxValueIndex = UINT16_MAX;

    bool trackingIncrease = false;         // Flag to indicate if we're currently tracking an increase
    int consecutiveDecreaseCount = 0;      // Counter for consecutive decreases within an increasing trend
    double currentMaxGradient = -INFINITY; // Maximum gradient within the current increasing trend
    uint16_t currentMaxGradientIndex = UINT16_MAX; // Index of the maximum gradient within the current increasing trend

    DEBUG_PRINT_1("Threshold for Increase: %.6f\n", threshold);

    // **Iterate through the gradient array within the specified window size**
    for (uint16_t i = 0; i < windowSize; ++i) {
        double gradient = gradients[i];    // Current gradient value
        uint16_t currentIndex = i;         // Current index in the gradient array

        // **Ignore NaN values to prevent invalid calculations**
        if (isnan(gradient)) continue;

        DEBUG_PRINT_2("Gradient at index %u: %.6f\n", currentIndex, gradient);

        // **Check if the current gradient exceeds the specified threshold**
        if (gradient > threshold) {  // Significant positive gradient
            if (!trackingIncrease) {
                // **Start tracking a new increasing trend**

                // Set the start and end indices to the current index
                currentIncreaseInfo.startIndex = currentIndex;
                currentIncreaseInfo.endIndex = currentIndex;

                trackingIncrease = true;               // Set the flag indicating we're tracking an increase
                consecutiveDecreaseCount = 0;          // Reset the consecutive decrease counter
                currentMaxGradient = gradient;         // Initialize the current maximum gradient
                currentMaxGradientIndex = currentIndex; // Initialize the index of the maximum gradient
                currentIncreaseInfo.sumGradients = gradient; // Initialize the sum of gradients

                DEBUG_PRINT_3("Started new increasing trend at index %u\n", currentIndex);
            } else {
                // **Continue tracking the current increasing trend**

                currentIncreaseInfo.endIndex = currentIndex; // Update the end index of the current trend

                // **Update the maximum gradient within the current trend**
                if (gradient > currentMaxGradient) {
                    currentMaxGradient = gradient;
                    currentMaxGradientIndex = currentIndex;
                }

                // **Reset the consecutive decrease counter since we have an increase**
                consecutiveDecreaseCount = 0;

                // **Add the current gradient to the sum of gradients**
                currentIncreaseInfo.sumGradients += gradient;
            }

        } else {
            // **Handle gradients that do not exceed the threshold**

            if (trackingIncrease) {
                // **We are currently tracking an increase, but encountered a decrease or insignificant gradient**

                consecutiveDecreaseCount++; // Increment the consecutive decrease counter

                // **Add the current gradient to the sum, even if it's a decrease**
                currentIncreaseInfo.sumGradients += gradient;

                // **Check if the consecutive decreases exceed the allowed limit**
                if (consecutiveDecreaseCount > MAX_ALLOWED_DECREASES) {
                    // **Consider ending the current increasing trend**

                    // **Calculate the length of the current increasing trend**
                    uint16_t currentTrendLength = currentIncreaseInfo.endIndex - currentIncreaseInfo.startIndex + 1;

                    // **Calculate the length of the longest increasing trend found so far**
                    uint16_t longestTrendLength = (longestIncreaseInfo.startIndex != UINT16_MAX) ?
                        (longestIncreaseInfo.endIndex - longestIncreaseInfo.startIndex + 1) : 0;

                    // **Compare the lengths to determine if the current trend is the longest**
                    if (currentTrendLength > longestTrendLength) {
                        // **Update the longest increasing trend information**
                        longestIncreaseInfo = currentIncreaseInfo;
                        longestIncreaseInfo.maxValue = currentMaxGradient;
                        longestIncreaseInfo.maxValueIndex = currentMaxGradientIndex;

                        DEBUG_PRINT_1("New longest increasing trend from index %u to %u\n",
                                     longestIncreaseInfo.startIndex, longestIncreaseInfo.endIndex);
                    }

                    // **Stop tracking the current increasing trend**
                    trackingIncrease = false;
                    consecutiveDecreaseCount = 0; // Reset the counter

                    // **Reset the current increasing trend information**
                    memset(&currentIncreaseInfo, 0, sizeof(GradientTrendIndices));
                    currentIncreaseInfo.startIndex = UINT16_MAX;
                    currentIncreaseInfo.endIndex = UINT16_MAX;
                    currentIncreaseInfo.maxValueIndex = UINT16_MAX;
                }
            }
            // **If not tracking an increase, do nothing for decreases**
        }
    }

    // **After iterating, check if the current trend is the longest**
    if (trackingIncrease) {
        // **Calculate the length of the current increasing trend**
        uint16_t currentTrendLength = currentIncreaseInfo.endIndex - currentIncreaseInfo.startIndex + 1;

        // **Calculate the length of the longest increasing trend found so far**
        uint16_t longestTrendLength = (longestIncreaseInfo.startIndex != UINT16_MAX) ?
            (longestIncreaseInfo.endIndex - longestIncreaseInfo.startIndex + 1) : 0;

        // **Compare the lengths to determine if the current trend is the longest**
        if (currentTrendLength > longestTrendLength) {
            // **Update the longest increasing trend information**
            longestIncreaseInfo = currentIncreaseInfo;
            longestIncreaseInfo.maxValue = currentMaxGradient;
            longestIncreaseInfo.maxValueIndex = currentMaxGradientIndex;

            DEBUG_PRINT_2("Final longest increasing trend from index %u to %u\n",
                   longestIncreaseInfo.startIndex, longestIncreaseInfo.endIndex);
        }
    }

    DEBUG_PRINT_3("Exiting findConsistentIncreaseInternal\n");

    // **Return the longest increasing trend information**
    return longestIncreaseInfo;
}


/**
 * @brief Finds the longest consistent decreasing trend in the gradient array based on a specified threshold.
 *
 * This function scans through an array of gradient values to identify the longest consistent decreasing trend
 * that exceeds a specified threshold. It allows a specified number of minor consecutive increases to account
 * for small fluctuations. The function returns a `GradientTrendIndices` struct containing the start and end
 * indices of the longest consistent decrease, the minimum gradient value within that trend,
 * and the sum of gradients within the trend.
 *
 * @param gradients The array of gradient values to analyze.
 * @param startIndex The starting index in the original data array.
 * @param windowSize The number of elements in the gradient array to analyze.
 * @param threshold The threshold below which gradients are considered significant.
 * @return GradientTrendIndices A structure containing the trend information.
 */
/**
 * @brief Finds the longest consistent decreasing trend in the gradient array based on a specified threshold.
 *
 * This function scans through an array of gradient values to identify the longest consistent decreasing trend
 * that exceeds a specified threshold. It allows a specified number of minor consecutive increases to account
 * for small fluctuations. The function returns a `GradientTrendIndices` struct containing the start and end
 * indices of the longest consistent decrease, the minimum gradient value within that trend,
 * and the sum of gradients within the trend.
 *
 * @param gradients   The array of gradient values to analyze.
 * @param windowSize  The number of elements in the gradient array to analyze.
 * @param threshold   The threshold below which gradients are considered significant.
 * @return GradientTrendIndices A structure containing the trend information.
 */
static GradientTrendIndices findConsistentDecreaseInternal(
    const double *gradients,
    uint16_t windowSize,
    double threshold
) {
    DEBUG_PRINT_3("Entering findConsistentDecreaseInternal\n");

    // **Initialize the structure to store the longest decreasing trend found**
    GradientTrendIndices longestDecreaseInfo;
    memset(&longestDecreaseInfo, 0, sizeof(GradientTrendIndices));
    longestDecreaseInfo.startIndex = UINT16_MAX;      // Initialize to UINT16_MAX to indicate invalid index
    longestDecreaseInfo.endIndex = UINT16_MAX;
    longestDecreaseInfo.maxValueIndex = UINT16_MAX;

    // **Initialize the structure to store the current decreasing trend being tracked**
    GradientTrendIndices currentDecreaseInfo;
    memset(&currentDecreaseInfo, 0, sizeof(GradientTrendIndices));
    currentDecreaseInfo.startIndex = UINT16_MAX;
    currentDecreaseInfo.endIndex = UINT16_MAX;
    currentDecreaseInfo.maxValueIndex = UINT16_MAX;

    bool trackingDecrease = false;          // Flag to indicate if we're currently tracking a decrease
    int consecutiveIncreaseCount = 0;       // Counter for consecutive increases within a decreasing trend
    double currentMinGradient = INFINITY;   // Minimum gradient within the current decreasing trend
    uint16_t currentMinGradientIndex = UINT16_MAX; // Index of the minimum gradient within the current decreasing trend

    DEBUG_PRINT_1("Threshold for Decrease: %.6f\n", threshold);

    // **Iterate through the gradient array within the specified window size**
    for (uint16_t i = 0; i < windowSize; ++i) {
        double gradient = gradients[i];    // Current gradient value
        uint16_t currentIndex = i;         // Current index in the gradient array

        // **Ignore NaN values to prevent invalid calculations**
        if (isnan(gradient)) continue;

        DEBUG_PRINT_2("Gradient at index %u: %.6f\n", currentIndex, gradient);

        // **Check if the current gradient is less than the specified threshold (negative)**
        if (gradient < threshold) {  // Significant negative gradient
            if (!trackingDecrease) {
                // **Start tracking a new decreasing trend**

                // Set the start and end indices to the current index
                currentDecreaseInfo.startIndex = currentIndex;
                currentDecreaseInfo.endIndex = currentIndex;

                trackingDecrease = true;                // Set the flag indicating we're tracking a decrease
                consecutiveIncreaseCount = 0;           // Reset the consecutive increase counter
                currentMinGradient = gradient;          // Initialize the current minimum gradient
                currentMinGradientIndex = currentIndex; // Initialize the index of the minimum gradient
                currentDecreaseInfo.sumGradients = gradient; // Initialize the sum of gradients

                DEBUG_PRINT_3("Started new decreasing trend at index %u\n", currentIndex);
            } else {
                // **Continue tracking the current decreasing trend**

                currentDecreaseInfo.endIndex = currentIndex; // Update the end index of the current trend

                // **Update the minimum gradient within the current trend**
                if (gradient < currentMinGradient) {
                    currentMinGradient = gradient;
                    currentMinGradientIndex = currentIndex;
                }

                // **Reset the consecutive increase counter since we have a decrease**
                consecutiveIncreaseCount = 0;

                // **Add the current gradient to the sum of gradients**
                currentDecreaseInfo.sumGradients += gradient;
            }

        } else {
            // **Handle gradients that do not exceed the negative threshold**

            if (trackingDecrease) {
                // **We are currently tracking a decrease, but encountered an increase or insignificant gradient**

                consecutiveIncreaseCount++; // Increment the consecutive increase counter

                // **Add the current gradient to the sum, even if it's an increase**
                currentDecreaseInfo.sumGradients += gradient;

                // **Check if the consecutive increases exceed the allowed limit**
                if (consecutiveIncreaseCount > MAX_ALLOWED_INCREASES) {
                    // **Consider ending the current decreasing trend**

                    // **Calculate the length of the current decreasing trend**
                    uint16_t currentTrendLength = currentDecreaseInfo.endIndex - currentDecreaseInfo.startIndex + 1;

                    // **Calculate the length of the longest decreasing trend found so far**
                    uint16_t longestTrendLength = (longestDecreaseInfo.startIndex != UINT16_MAX) ?
                        (longestDecreaseInfo.endIndex - longestDecreaseInfo.startIndex + 1) : 0;

                    // **Compare the lengths to determine if the current trend is the longest**
                    if (currentTrendLength > longestTrendLength) {
                        // **Update the longest decreasing trend information**
                        longestDecreaseInfo = currentDecreaseInfo;
                        longestDecreaseInfo.maxValue = currentMinGradient;
                        longestDecreaseInfo.maxValueIndex = currentMinGradientIndex;

                        DEBUG_PRINT_1("New longest decreasing trend from index %u to %u\n",
                                     longestDecreaseInfo.startIndex, longestDecreaseInfo.endIndex);
                    }

                    // **Stop tracking the current decreasing trend**
                    trackingDecrease = false;
                    consecutiveIncreaseCount = 0; // Reset the counter

                    // **Reset the current decreasing trend information**
                    memset(&currentDecreaseInfo, 0, sizeof(GradientTrendIndices));
                    currentDecreaseInfo.startIndex = UINT16_MAX;
                    currentDecreaseInfo.endIndex = UINT16_MAX;
                    currentDecreaseInfo.maxValueIndex = UINT16_MAX;
                }
            }
            // **If not tracking a decrease, do nothing for increases**
        }
    }

    // **After iterating, check if the current trend is the longest**
    if (trackingDecrease) {
        // **Calculate the length of the current decreasing trend**
        uint16_t currentTrendLength = currentDecreaseInfo.endIndex - currentDecreaseInfo.startIndex + 1;

        // **Calculate the length of the longest decreasing trend found so far**
        uint16_t longestTrendLength = (longestDecreaseInfo.startIndex != UINT16_MAX) ?
            (longestDecreaseInfo.endIndex - longestDecreaseInfo.startIndex + 1) : 0;

        // **Compare the lengths to determine if the current trend is the longest**
        if (currentTrendLength > longestTrendLength) {
            // **Update the longest decreasing trend information**
            longestDecreaseInfo = currentDecreaseInfo;
            longestDecreaseInfo.maxValue = currentMinGradient;
            longestDecreaseInfo.maxValueIndex = currentMinGradientIndex;

            DEBUG_PRINT_2("Final longest decreasing trend from index %u to %u\n",
                   longestDecreaseInfo.startIndex, longestDecreaseInfo.endIndex);
        }
    }

    DEBUG_PRINT_3("Exiting findConsistentDecreaseInternal\n");

    // **Return the longest decreasing trend information**
    return longestDecreaseInfo;
}


/**
 * @brief Determines the move direction based on consistent increase and decrease counts.
 *
 * This function evaluates the consistent increase and decrease counts from the gradient trends
 * and the gradient array to determine the move direction in the sliding window analysis.
 * It follows a series of logical steps, using thresholds and additional checks, to decide
 * whether to move left, right, or consider the position as on the peak or undecided.
 *
 * @param longestIncrease Pointer to the longest consistent increasing trend indices.
 * @param longestDecrease Pointer to the longest consistent decreasing trend indices.
 * @param gradients       The array of gradient values to analyze.
 * @param gradientSize    The number of elements in the gradient array.
 * @return PeakPosition The determined move direction.
 */
PeakPosition determineMoveDirection(
    const GradientTrendIndices* longestIncrease,
    const GradientTrendIndices* longestDecrease,
    const double* gradients,
    uint16_t gradientSize
) {
    // **Adjust gradientSize and gradients to ignore the first value**
    const double* adjustedGradients = &gradients[1];
    uint16_t adjustedGradientSize = gradientSize - 1;

    // **Step 1: Calculate the counts of consistent increases and decreases**
    // Initialize increaseCount and decreaseCount to zero.
    uint16_t increaseCount = 0;
    uint16_t decreaseCount = 0;

    // Adjust the longest increase trend indices if they include the first gradient value
    uint16_t increaseStartIndex = longestIncrease->startIndex;
    uint16_t increaseEndIndex = longestIncrease->endIndex;

    if (increaseStartIndex == 0) {
        increaseStartIndex = 1;
    }

    // Calculate the number of points in the consistent increasing trend
    if (increaseStartIndex != UINT16_MAX && increaseEndIndex != UINT16_MAX && increaseEndIndex >= increaseStartIndex) {
        increaseCount = increaseEndIndex - increaseStartIndex + 1;
    }

    // Adjust the longest decrease trend indices if they include the first gradient value
    uint16_t decreaseStartIndex = longestDecrease->startIndex;
    uint16_t decreaseEndIndex = longestDecrease->endIndex;

    if (decreaseStartIndex == 0) {
        decreaseStartIndex = 1;
    }

    // Calculate the number of points in the consistent decreasing trend
    if (decreaseStartIndex != UINT16_MAX && decreaseEndIndex != UINT16_MAX && decreaseEndIndex >= decreaseStartIndex) {
        decreaseCount = decreaseEndIndex - decreaseStartIndex + 1;
    }

    // **Debugging Output**
    printf("Adjusted Increase Count: %u, Adjusted Decrease Count: %u\n", increaseCount, decreaseCount);
    //DEBUG_PRINT_2("Adjusted Longest Increase: startIndex=%u, endIndex=%u\n", increaseStartIndex, increaseEndIndex);
    //DEBUG_PRINT_2("Adjusted Longest Decrease: startIndex=%u, endIndex=%u\n", decreaseStartIndex, decreaseEndIndex);

    // **Step 2: Define the trend threshold**
    uint16_t trendThreshold = TREND_THRESHOLD; // e.g., TREND_THRESHOLD = RLS_WINDOW / 4

    // **Debugging Output**
    DEBUG_PRINT_1("Trend Threshold: %u\n", trendThreshold);

    // **Step 3: Initialize the moveDirection to UNDECIDED**
    PeakPosition moveDirection = UNDECIDED;

    // **Step 4: Analyze the maximum and minimum gradient values**
    double maxGradient = -INFINITY;
    double minGradient = INFINITY;

    // Iterate over the adjusted gradient array to find the global max and min gradients
    for (uint16_t i = 0; i < adjustedGradientSize; ++i) {
        double gradient = adjustedGradients[i];
        if (gradient > maxGradient) {
            maxGradient = gradient;
        }
        if (gradient < minGradient) {
            minGradient = gradient;
        }
    }

    // **Debugging Output**
    DEBUG_PRINT_1("Max Gradient: %.6f, Min Gradient: %.6f\n", maxGradient, minGradient);

    // **Step 5: Determine if gradients are within undecided thresholds**
    bool gradientsWithinUndecidedThresholds = (maxGradient < MAX_GRADIENT_THRESHOLD && minGradient > MIN_GRADIENT_THRESHOLD);

    // **Step 6: Decide moveDirection based on counts and gradients**
    if (!gradientsWithinUndecidedThresholds && (increaseCount >= trendThreshold || decreaseCount >= trendThreshold)) {
        // **Sub-case 1: Significant increase count and gradients**
        if (increaseCount >= trendThreshold && decreaseCount < trendThreshold) {
            printf("Condition Met: Significant increase count and gradients.\n");
            moveDirection = RIGHT_SIDE;
        }
        // **Sub-case 2: Significant decrease count and gradients**
        else if (decreaseCount >= trendThreshold && increaseCount < trendThreshold) {
            printf("Condition Met: Significant decrease count and gradients.\n");
            moveDirection = LEFT_SIDE;
        }
        // **Sub-case 3: Both counts are significant**
        else if (increaseCount >= trendThreshold && decreaseCount >= trendThreshold) {                         // TAMAMEN SAÇMA SAPAN BİR LOGIC. 
           
           // Compare the sum of gradients from the passed-in structs
            double sumIncreaseGradients = longestIncrease->sumGradients;
            double sumDecreaseGradients = longestDecrease->sumGradients;

            if (fabs(sumDecreaseGradients) > fabs(sumIncreaseGradients)) {
                printf("Decrease trend sum is higher. Moving LEFT.\n");
                moveDirection = LEFT_SIDE;
            } else {
                printf("Increase trend sum is higher. Moving RIGHT.\n");
                moveDirection = RIGHT_SIDE;
            }
        }
    } else {
        DEBUG_PRINT_1("Counts below threshold or gradients within undecided thresholds; proceeding to additional checks.\n");
        
        //thresholdu açan

        // **Sub-step 6.1: Check if gradients are within undecided thresholds**
        if (gradientsWithinUndecidedThresholds) {
            printf("All gradients are within undecided thresholds (%.2f, %.2f).\n", MIN_GRADIENT_THRESHOLD, MAX_GRADIENT_THRESHOLD);
            moveDirection = UNDECIDED;
        }
        // **Sub-step 6.2: Compare gradients in left and right halves**
        else {
            DEBUG_PRINT_2("Gradients exceed undecided thresholds; comparing halves.\n");

            // Split the adjusted gradient array into two halves
            uint16_t midIndex = adjustedGradientSize / 2;
            double leftMax = -INFINITY;
            double rightMax = -INFINITY;

            // Find the maximum gradient in the left half
            for (uint16_t i = 0; i < midIndex; ++i) {
                if (adjustedGradients[i] > leftMax) {
                    leftMax = adjustedGradients[i];
                }
            }

            // Find the maximum gradient in the right half
            for (uint16_t i = midIndex; i < adjustedGradientSize; ++i) {
                if (adjustedGradients[i] > rightMax) {
                    rightMax = adjustedGradients[i];
                }
            }

            // **Debugging Output**
            DEBUG_PRINT_2("Left Max Gradient: %.6f, Right Max Gradient: %.6f\n", leftMax, rightMax);

            // **Sub-step 6.3: Decide moveDirection based on the comparison**
            if (rightMax > leftMax) {
                DEBUG_PRINT_2("Right half has higher max gradient.\n");
                moveDirection = RIGHT_SIDE;
            } else if (leftMax > rightMax) {
                DEBUG_PRINT_2("Left half has higher max gradient.\n");
                moveDirection = LEFT_SIDE;
            } else {
                DEBUG_PRINT_2("Both halves have equal max gradients.\n");
                moveDirection = UNDECIDED;
            }
        }
    }
    // **Step 9: Return the determined moveDirection**
    return moveDirection;
}


/**
 * @brief Identifies gradient trends within a window using RLS regression and MAD-based statistics.
 *
 * This function calculates gradients using the provided `calculate_gradient` function, identifies the
 * longest consistent increasing and decreasing trends using `findConsistentIncreaseInternal` and
 * `findConsistentDecreaseInternal`, computes separate statistical measures (median and MAD)
 * for each trend, and stores the results in the provided `gradCalcResultInc` and `gradCalcResultDec` structures.
 *
 * @param data              Array of data points (MqsRawDataPoint_t) to analyze.
 * @param startIndex        The starting index in the data array.
 * @param analysisLength    The number of data points to include in the gradient analysis window.
 * @param degree            The degree of the polynomial used for regression.
 * @param calculate_gradient Function pointer to the gradient calculation function.
 * @param gradCalcResultInc Pointer to store the gradient calculation results for increases.
 * @param gradCalcResultDec Pointer to store the gradient calculation results for decreases.
 *
 * @return PeakPosition The move direction (LEFT_SIDE, RIGHT_SIDE, ON_PEAK, or UNDECIDED).
 */
PeakPosition identifyTrends(
    const MqsRawDataPoint_t *data,
    uint16_t startIndex,
    uint16_t analysisLength,
    uint8_t degree,
    double (*calculate_gradient)(const RegressionState *, double),
    GradientCalculationResult *gradCalcResultInc,
    GradientCalculationResult *gradCalcResultDec
) {
    // Initialize gradCalcResult structures
    memset(gradCalcResultInc, 0, sizeof(GradientCalculationResult));
    memset(gradCalcResultDec, 0, sizeof(GradientCalculationResult));

    DEBUG_PRINT_1("Starting identifyTrends with startIndex=%u, analysisLength=%u\n",
                  startIndex, analysisLength);

    // Calculate all gradients over this window
    GradientCalculationResult gradCalcResultAll;
    memset(&gradCalcResultAll, 0, sizeof(GradientCalculationResult));

    // Debugging: Calling trackGradients to calculate gradients
    DEBUG_PRINT_2("Calling trackGradients (MqsRawDataPoint_t) to calculate gradients.\n");
    trackGradients(
        data,               // Now MqsRawDataPoint_t array
        analysisLength,
        startIndex,
        degree,
        calculate_gradient,
        &gradCalcResultAll
    );

    uint16_t count = gradCalcResultAll.size;
    DEBUG_PRINT_1("Total gradients calculated: %u\n", count);

    if (count == 0) {
        DEBUG_PRINT_1("No gradients calculated. Exiting identifyTrends.\n");
        return UNDECIDED;
    }

    // Find the longest consistent increasing trend
    DEBUG_PRINT_1("Finding the longest consistent increasing trend.\n");
    GradientTrendIndices longestIncrease = findConsistentIncreaseInternal(
        gradCalcResultAll.gradients, // array of double gradients
        count,
        0.0 // Adjust threshold as needed
    );

    // Adjust indices to absolute positions
    if (longestIncrease.startIndex != UINT16_MAX) {
        DEBUG_PRINT_2("Adjusting indices of the longest increasing trend to absolute.\n");
        DEBUG_PRINT_3("Before: startIndex=%u, endIndex=%u, maxValueIndex=%u\n",
                      longestIncrease.startIndex, longestIncrease.endIndex, longestIncrease.maxValueIndex);

        longestIncrease.startIndex   += startIndex;
        longestIncrease.endIndex     += startIndex;
        longestIncrease.maxValueIndex+= startIndex;

        DEBUG_PRINT_3("After:  startIndex=%u, endIndex=%u, maxValueIndex=%u\n",
                      longestIncrease.startIndex, longestIncrease.endIndex, longestIncrease.maxValueIndex);
    } else {
        DEBUG_PRINT_1("No increasing trend found.\n");
    }

    // Calculate median & MAD for the longest increasing trend
    if (longestIncrease.startIndex != UINT16_MAX && longestIncrease.endIndex != UINT16_MAX) {
        uint16_t relativeStart = longestIncrease.startIndex - startIndex;
        uint16_t incLength     = longestIncrease.endIndex - longestIncrease.startIndex + 1;

        DEBUG_PRINT_2("Calculating median & MAD for the longest increasing trend.\n");
        DEBUG_PRINT_3("Relative start index: %u\n", relativeStart);
        DEBUG_PRINT_3("Length of inc trend:  %u\n", incLength);

        if (incLength > 0 && (relativeStart + incLength <= count)) {
            gradCalcResultInc->size = incLength;
            memcpy(gradCalcResultInc->gradients,
                   &gradCalcResultAll.gradients[relativeStart],
                   incLength * sizeof(double));

            calculateMedianAndMAD(
                gradCalcResultInc->gradients,
                gradCalcResultInc->size,
                &gradCalcResultInc->median,
                &gradCalcResultInc->mad
            );

            gradCalcResultInc->startIndex = longestIncrease.startIndex;
            gradCalcResultInc->endIndex   = longestIncrease.endIndex;

            DEBUG_PRINT_2("Median (inc) = %.6f, MAD (inc) = %.6f\n",
                          gradCalcResultInc->median,
                          gradCalcResultInc->mad);
        } else {
            DEBUG_PRINT_1("Invalid inc indices. Skipping median & MAD.\n");
        }
    }

    // Find the longest consistent decreasing trend
    DEBUG_PRINT_1("Finding the longest consistent decreasing trend.\n");
    GradientTrendIndices longestDecrease = findConsistentDecreaseInternal(
        gradCalcResultAll.gradients,
        count,
        0.0 // Adjust threshold as needed
    );

    // Adjust indices to absolute positions
    if (longestDecrease.startIndex != UINT16_MAX) {
        DEBUG_PRINT_2("Adjusting indices of the longest decreasing trend to absolute.\n");
        DEBUG_PRINT_3("Before: startIndex=%u, endIndex=%u, maxValueIndex=%u\n",
                      longestDecrease.startIndex, longestDecrease.endIndex, longestDecrease.maxValueIndex);

        longestDecrease.startIndex    += startIndex;
        longestDecrease.endIndex      += startIndex;
        longestDecrease.maxValueIndex += startIndex;

        DEBUG_PRINT_3("After:  startIndex=%u, endIndex=%u, maxValueIndex=%u\n",
                      longestDecrease.startIndex, longestDecrease.endIndex, longestDecrease.maxValueIndex);
    } else {
        DEBUG_PRINT_1("No decreasing trend found.\n");
    }

    // Calculate median & MAD for the longest decreasing trend
    if (longestDecrease.startIndex != UINT16_MAX && longestDecrease.endIndex != UINT16_MAX) {
        uint16_t relativeStart = longestDecrease.startIndex - startIndex;
        uint16_t decLength     = longestDecrease.endIndex - longestDecrease.startIndex + 1;

        DEBUG_PRINT_2("Calculating median & MAD for the longest decreasing trend.\n");
        DEBUG_PRINT_3("Relative start index: %u\n", relativeStart);
        DEBUG_PRINT_3("Length of dec trend:  %u\n", decLength);

        if (decLength > 0 && (relativeStart + decLength <= count)) {
            gradCalcResultDec->size = decLength;
            memcpy(gradCalcResultDec->gradients,
                   &gradCalcResultAll.gradients[relativeStart],
                   decLength * sizeof(double));

            calculateMedianAndMAD(
                gradCalcResultDec->gradients,
                gradCalcResultDec->size,
                &gradCalcResultDec->median,
                &gradCalcResultDec->mad
            );

            gradCalcResultDec->startIndex = longestDecrease.startIndex;
            gradCalcResultDec->endIndex   = longestDecrease.endIndex;

            DEBUG_PRINT_2("Median (dec) = %.6f, MAD (dec) = %.6f\n",
                          gradCalcResultDec->median,
                          gradCalcResultDec->mad);
        } else {
            DEBUG_PRINT_1("Invalid dec indices. Skipping median & MAD.\n");
        }
    }
    
    // Finally determine move direction
    PeakPosition moveDirection = determineMoveDirection(
        &longestIncrease,
        &longestDecrease,
        gradCalcResultAll.gradients, // the array of computed gradients
        count
    );
    
    // Optionally, log or store the moveDirection
    
    DEBUG_PRINT_1("Completed identifyTrends.\n");
    return moveDirection;
}

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
) {
    int32_t newStartIndex = (int32_t)currentStartIndex + moveAmount;

    if (newStartIndex < 0) {
        newStartIndex = 0;
    } else if ((uint32_t)newStartIndex + windowSize > dataLength) {
        if (dataLength >= windowSize) {
            newStartIndex = dataLength - windowSize;
        } else {
            newStartIndex = 0;  // If dataLength < windowSize, start at 0
        }
    }

    return (uint16_t)newStartIndex;
}

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
) {
    DEBUG_PRINT_3("Detecting trends.\n");
    
    // Initialize thresholds
    double thresholdIncrease = 0.0;
    double thresholdDecrease = 0.0;
    
    // Depending on trendType, set thresholds
    if (trendType == TREND_TYPE_ABSOLUTE) {
        thresholdIncrease = gradCalcResultInc->median;
        thresholdDecrease = gradCalcResultDec->median;
    } else if (trendType == TREND_TYPE_SIGNIFICANT) {
        thresholdIncrease = gradCalcResultInc->median + (2 * gradCalcResultInc->mad);
        thresholdDecrease = gradCalcResultDec->median - (2 * gradCalcResultDec->mad);
    } else {
        DEBUG_PRINT_1("Invalid trend detection type specified.\n");
        return;
    }
    
    // Initialize move amount
    *(trendResultBase->moveAmount) = 0;
    
    // Detect increases
    if (gradCalcResultInc->size > 0) {
        for (uint16_t i = 0; i < gradCalcResultInc->size; ++i) {
            double gradient = gradCalcResultInc->gradients[i];
            if (gradient >= thresholdIncrease) {
                // Detected an increase
                if (trendResultBase->increase->startIndex == UINT16_MAX &&
                    trendResultBase->increase->endIndex == UINT16_MAX &&
                    trendResultBase->increase->sumGradients == 0.0) {
                    trendResultBase->increase->startIndex = i;
                    trendResultBase->increase->endIndex = i;
                    trendResultBase->increase->maxValue = gradient;
                    trendResultBase->increase->maxValueIndex = i;
                    trendResultBase->increase->sumGradients = gradient;
                } else {
                    trendResultBase->increase->endIndex = i;
                    trendResultBase->increase->sumGradients += gradient;
                    if (gradient > trendResultBase->increase->maxValue) {
                        trendResultBase->increase->maxValue = gradient;
                        trendResultBase->increase->maxValueIndex = i;
                    }
                }
    
                // Increment move amount
                (*(trendResultBase->moveAmount))++;
    
                DEBUG_PRINT_3("Detected increase at index %u.\n", i);
            }
        }
    }
    
    // Detect decreases
    if (gradCalcResultDec->size > 0) {
        for (uint16_t i = 0; i < gradCalcResultDec->size; ++i) {
            double gradient = gradCalcResultDec->gradients[i];
            if (gradient <= thresholdDecrease) {
                // Detected a decrease
                if (trendResultBase->decrease->startIndex == UINT16_MAX &&
                    trendResultBase->decrease->endIndex == UINT16_MAX &&
                    trendResultBase->decrease->sumGradients == 0.0) {
                    trendResultBase->decrease->startIndex = i;
                    trendResultBase->decrease->endIndex = i;
                    trendResultBase->decrease->maxValue = gradient;
                    trendResultBase->decrease->maxValueIndex = i;
                    trendResultBase->decrease->sumGradients = gradient;
                } else {
                    trendResultBase->decrease->endIndex = i;
                    trendResultBase->decrease->sumGradients += gradient;
                    if (gradient < trendResultBase->decrease->maxValue) {
                        trendResultBase->decrease->maxValue = gradient;
                        trendResultBase->decrease->maxValueIndex = i;
                    }
                }
    
                // Increment move amount
                (*(trendResultBase->moveAmount))++;
    
                DEBUG_PRINT_3("Detected decrease at index %u.\n", i);
            }
        }
    }
    
    
    DEBUG_PRINT_2("Trend detection completed. Move amount: %d\n", *(trendResultBase->moveAmount));
}

/**
 * @brief Identifies gradient trends within a specified window.
 *
 * This function detects consistent increasing and decreasing trends based on first-order or second-order gradients
 * calculated using RLS regression. It delegates the detection to specific functions based on the
 * trend detection type (absolute, significant, or gradient trends).
 *
 * @param data             Array of data points (MqsRawDataPoint_t) to analyze.
 * @param startIndex       The starting index in the data array.
 * @param analysisLength   The number of data points to include in the gradient analysis window.
 * @param degree           The degree of the polynomial used for regression.
 * @param gradientOrder    The order of the gradient (first-order or second-order).
 * @param trendType        The type of trend detection to perform.
 * @param trendResult      Pointer to store the trend detection result.
 * @param gradCalcResultInc Pointer to store the gradient calculation results for increases.
 * @param gradCalcResultDec Pointer to store the gradient calculation results for decreases.
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
) {
    // Determine the gradient calculation function based on gradientOrder
    double (*calculate_gradient)(const RegressionState *, double) = NULL;
// Initialize only if debug level is greater than 0
#if DEBUG_LEVEL > 0
    const char* gradientOrderStr = NULL;
#endif

switch (gradientOrder) {
    case GRADIENT_ORDER_FIRST:
        calculate_gradient = calculate_first_order_gradient;
#if DEBUG_LEVEL > 0
        gradientOrderStr   = "First-Order";
#endif
        break;
    case GRADIENT_ORDER_SECOND:
        calculate_gradient = calculate_second_order_gradient;
#if DEBUG_LEVEL > 0
        gradientOrderStr   = "Second-Order";
#endif
        break;
    default:
        DEBUG_PRINT_1("Invalid gradient order specified.\n");
        return;
}

    // Print or log which gradient order
    DEBUG_PRINT_2("=== Calculating %s Gradients ===\n", gradientOrderStr);

    // We call identifyTrends with MqsRawDataPoint_t
    PeakPosition moveDirection = identifyTrends(
        data,
        startIndex,
        analysisLength,
        degree,
        calculate_gradient,
        gradCalcResultInc,
        gradCalcResultDec
    );

    // If no gradients, skip further processing
    if (gradCalcResultInc->size == 0 && gradCalcResultDec->size == 0) {
        DEBUG_PRINT_1("No gradients detected, skipping trend detection.\n");
        return;
    }

    // Handle different trend detection types
    switch (trendType) {
        case TREND_TYPE_ABSOLUTE:
        case TREND_TYPE_SIGNIFICANT: {
            GradientTrendResultBase *trendResultBase = NULL;

            if (trendType == TREND_TYPE_ABSOLUTE) {
                GradientTrendResultAbsolute *absResult = (GradientTrendResultAbsolute *)trendResult;
                initializeGradientTrendResultAbsolute(absResult);
                trendResultBase = &(absResult->base);
                DEBUG_PRINT_1("=== Detecting Absolute %s Gradient Trends ===\n", gradientOrderStr);
            } else { // TREND_TYPE_SIGNIFICANT
                GradientTrendResultSignificant *sigResult = (GradientTrendResultSignificant *)trendResult;
                initializeGradientTrendResultSignificant(sigResult);
                trendResultBase = &(sigResult->base);
                DEBUG_PRINT_1("=== Detecting Significant %s Gradient Trends ===\n", gradientOrderStr);
            }

            // We already have moveDirection from identifyTrends
            *(trendResultBase->moveDirection) = moveDirection;

            // Print results depending on the type
            if (trendType == TREND_TYPE_ABSOLUTE) {
                GradientTrendResultAbsolute *absResult = (GradientTrendResultAbsolute *)trendResult;
                DEBUG_PRINT_1("Absolute %s Gradient Trends:\n", gradientOrderStr);
                DEBUG_PRINT_1("Move amount: %d\n", absResult->moveAmountAbsolute);
                printTrendDetails("Absolute Increase", &absResult->absoluteIncrease);
                printTrendDetails("Absolute Decrease", &absResult->absoluteDecrease);
            } else { // TREND_TYPE_SIGNIFICANT
                GradientTrendResultSignificant *sigResult = (GradientTrendResultSignificant *)trendResult;
                DEBUG_PRINT_1("Significant %s Gradient Trends:\n", gradientOrderStr);
                DEBUG_PRINT_1("Move amount: %d\n", sigResult->moveAmountSignificant);
                printTrendDetails("Significant Increase", &sigResult->significantIncrease);
                printTrendDetails("Significant Decrease", &sigResult->significantDecrease);
            }
            break;
        }

        case TREND_TYPE_GRADIENT_TRENDS: {
            GradientTrendResultAbsolute *gradResult = (GradientTrendResultAbsolute *)trendResult;
            initializeGradientTrendResultAbsolute(gradResult);

            DEBUG_PRINT_1("=== Detecting %s Gradient Trends ===\n", gradientOrderStr);
            DEBUG_PRINT_1("Step: Extracting consistent increases and decreases for %s gradients.\n", gradientOrderStr);

            // Build absoluteIncrease from gradCalcResultInc
            if (gradCalcResultInc->size > 0) {
                gradResult->absoluteIncrease.startIndex = gradCalcResultInc->startIndex;
                gradResult->absoluteIncrease.endIndex   = gradCalcResultInc->endIndex;

                // Find maxValue in the inc array
                double maxValue = -INFINITY;
                uint16_t maxValueIndex = UINT16_MAX;
                for (uint16_t i = 0; i < gradCalcResultInc->size; ++i) {
                    double g = gradCalcResultInc->gradients[i];
                    if (g > maxValue) {
                        maxValue       = g;
                        maxValueIndex  = gradCalcResultInc->startIndex + i;
                    }
                }
                gradResult->absoluteIncrease.maxValue      = maxValue;
                gradResult->absoluteIncrease.maxValueIndex = maxValueIndex;

                // Sum of inc gradients
                double sumGradients = 0.0;
                for (uint16_t i = 0; i < gradCalcResultInc->size; ++i) {
                    sumGradients += gradCalcResultInc->gradients[i];
                }
                gradResult->absoluteIncrease.sumGradients = sumGradients;
            } else {
                gradResult->absoluteIncrease.startIndex = UINT16_MAX;
                gradResult->absoluteIncrease.endIndex   = UINT16_MAX;
                gradResult->absoluteIncrease.maxValue   = 0.0;
                gradResult->absoluteIncrease.maxValueIndex = UINT16_MAX;
                gradResult->absoluteIncrease.sumGradients  = 0.0;
            }

            // Build absoluteDecrease from gradCalcResultDec
            if (gradCalcResultDec->size > 0) {
                gradResult->absoluteDecrease.startIndex = gradCalcResultDec->startIndex;
                gradResult->absoluteDecrease.endIndex   = gradCalcResultDec->endIndex;

                // Find minValue in the dec array
                double minValue = INFINITY;
                uint16_t minValueIndex = UINT16_MAX;
                for (uint16_t i = 0; i < gradCalcResultDec->size; ++i) {
                    double g = gradCalcResultDec->gradients[i];
                    if (g < minValue) {
                        minValue       = g;
                        minValueIndex  = gradCalcResultDec->startIndex + i;
                    }
                }
                gradResult->absoluteDecrease.maxValue      = minValue; // naming note: `maxValue` but used for min
                gradResult->absoluteDecrease.maxValueIndex = minValueIndex;

                // Sum of dec gradients
                double sumGradients = 0.0;
                for (uint16_t i = 0; i < gradCalcResultDec->size; ++i) {
                    sumGradients += gradCalcResultDec->gradients[i];
                }
                gradResult->absoluteDecrease.sumGradients = sumGradients;
            } else {
                gradResult->absoluteDecrease.startIndex = UINT16_MAX;
                gradResult->absoluteDecrease.endIndex   = UINT16_MAX;
                gradResult->absoluteDecrease.maxValue   = 0.0;
                gradResult->absoluteDecrease.maxValueIndex = UINT16_MAX;
                gradResult->absoluteDecrease.sumGradients  = 0.0;
            }

            // Now set the final moveDirection (already obtained from identifyTrends)
            gradResult->moveDirection     = moveDirection;
            gradResult->isSignificantPeak = (moveDirection == ON_PEAK);

            // Print or log
            printTrendDetails("Increase", &gradResult->absoluteIncrease);
            printTrendDetails("Decrease", &gradResult->absoluteDecrease);
            break;
        }

        default:
            DEBUG_PRINT_1("Invalid trend detection type specified.\n");
            break;
    }
}

// ** OBSOLETE FUNCTION ** 

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
) {
    // Initialize gradient calculation results
    GradientCalculationResult gradCalcResultAll;
    memset(&gradCalcResultAll, 0, sizeof(GradientCalculationResult));

    // Determine the gradient calculation function based on gradientOrder
    double (*calculate_gradient)(const RegressionState *, double) = NULL;
    const char* gradientOrderStr = NULL;

    switch (gradientOrder) {
        case GRADIENT_ORDER_FIRST:
            calculate_gradient = calculate_first_order_gradient;
            gradientOrderStr = "First-Order";
            break;
        case GRADIENT_ORDER_SECOND:
            calculate_gradient = calculate_second_order_gradient;
            gradientOrderStr = "Second-Order";
            break;
        default:
            DEBUG_PRINT_1("Invalid gradient order specified.\n");
            return NO_TREND;
    }

    printf("=== Calculating %s Gradients ===\n", gradientOrderStr);

    // Calculate gradients using trackGradients
    trackGradients(
        data,
        analysisLength,
        startIndex,
        degree,
        calculate_gradient,
        &gradCalcResultAll
    );

    // Check if any gradients have been calculated
    if (gradCalcResultAll.size == 0) {
        DEBUG_PRINT_1("No gradients calculated. Exiting identifySplitMoveDirection.\n");
        return NO_TREND;
    }

    printf("Total gradients calculated: %u\n", gradCalcResultAll.size);

    // Split the gradient array into two halves
    uint16_t midIndex = gradCalcResultAll.size / 2;
    uint16_t leftSize = midIndex;
    //uint16_t rightSize = gradCalcResultAll.size - midIndex;

    double leftSum = 0.0;
    double rightSum = 0.0;

    // Calculate sum of gradients in the left half
    for (uint16_t i = 0; i < leftSize; ++i) {
        leftSum += gradCalcResultAll.gradients[i];
    }

    // Calculate sum of gradients in the right half
    for (uint16_t i = midIndex; i < gradCalcResultAll.size; ++i) {
        rightSum += gradCalcResultAll.gradients[i];
    }

    DEBUG_PRINT_2("Left Sum: %.6f, Right Sum: %.6f\n", leftSum, rightSum);

    // Determine move direction based on sums
    if (leftSum > rightSum) {
        // Gradients are increasing towards the right; move window to the right
        printf("Determined Move Direction: MOVE_RIGHT\n");
        return MOVE_RIGHT;
    } else if (rightSum > leftSum) {
        // Gradients are increasing towards the left; move window to the left
        printf("Determined Move Direction: MOVE_LEFT\n");
        return MOVE_LEFT;
    } else {
        // Gradients are balanced; assume window is at the peak
        printf("Determined Move Direction: ON_THE_PEAK\n");
        return ON_THE_PEAK;
    }
}

// only 3 functions need MQsRawDataPoints changed. 

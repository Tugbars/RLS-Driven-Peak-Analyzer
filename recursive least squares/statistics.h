/**
 * @file statistics.h
 * @brief Header file for statistical functions.
 *
 * This module provides robust statistical functions such as median, percentile, and
 * median absolute deviation (MAD), which are essential for peak detection algorithms.
 * The functions are designed to handle outliers and provide reliable metrics for trend analysis.
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdint.h>
#include <stdbool.h>

/* Define maximum window size for analysis */
#define MAX_WINDOW_SIZE 30  // Adjust as necessary based on expected data sizes

/**
 * @brief Comparator function for qsort.
 *
 * Used internally for sorting double arrays.
 *
 * @param a Pointer to the first element.
 * @param b Pointer to the second element.
 * @return int Negative if *a < *b, positive if *a > *b, zero if equal.
 */
int compare_doubles(const void *a, const void *b);

/**
 * @brief Calculates the specified percentile of an array.
 *
 * This function sorts the data and computes the desired percentile value.
 *
 * @param data Array of double values. Note: The array will be modified (sorted).
 * @param size Number of elements in the array.
 * @param percentile Percentile to calculate (0-100).
 * @return double The calculated percentile value. Returns NAN if input is invalid.
 */
double calculate_percentile(double *data, uint16_t size, double percentile);

/**
 * @brief Calculates the median of an array.
 *
 * @param data Array of double values. Note: The array will be modified (sorted).
 * @param size Number of elements in the array.
 * @return double The median value. Returns NAN if input is invalid.
 */
double calculate_median(double *data, uint16_t size);

/**
 * @brief Calculates the Median Absolute Deviation (MAD) of an array.
 *
 * @param data Array of double values.
 * @param size Number of elements in the array.
 * @param median The median of the data array.
 * @return double The MAD value. Returns NAN if input is invalid.
 */
double calculate_mad(const double *data, uint16_t size, double median);

/**
 * @brief Calculates the median and MAD of gradient differences, excluding the top 1% percentile.
 *
 * This function filters out the top 1% of the data to eliminate potential noise peaks,
 * then calculates the median and MAD of the remaining data.
 *
 * @param data Array of gradient differences.
 * @param size Number of elements in the data array.
 * @param median Pointer to store the calculated median.
 * @param mad Pointer to store the calculated MAD.
 * @return bool Returns true if calculation was successful, false otherwise.
 */
bool calculateMedianAndMAD(const double *data, uint16_t size, double *median, double *mad);

#endif // STATISTICS_H

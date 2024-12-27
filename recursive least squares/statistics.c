/**
 * @file statistics.c
 * @brief Implementation of statistical functions.
 *
 * This module provides implementations for median, percentile, and MAD calculations,
 * which are essential for robust peak detection algorithms.
 */

#include "statistics.h"
#include <stdlib.h> // For qsort
#include <math.h>   // For NAN, fabs
#include <string.h> // For memcpy


/**
 * @brief Comparator function for qsort.
 */
int compare_doubles(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    else if (da > db) return 1;
    else return 0;
}

/**
 * @brief Calculates the specified percentile of an array.
 *
 * @param data Array of double values.
 * @param size Number of elements in the array.
 * @param percentile Percentile to calculate (0-100).
 * @return double The calculated percentile value.
 */
double calculate_percentile(double *data, uint16_t size, double percentile) {
    if (size == 0 || percentile < 0.0 || percentile > 100.0) return NAN;
    
    // Create a copy to avoid modifying the original data
    double data_copy[MAX_WINDOW_SIZE];
    if (size > MAX_WINDOW_SIZE) size = MAX_WINDOW_SIZE; // Prevent buffer overflow
    memcpy(data_copy, data, sizeof(double) * size);
    
    qsort(data_copy, size, sizeof(double), compare_doubles);
    
    double pos = (percentile / 100.0) * (size - 1);
    uint16_t index = (uint16_t)floor(pos);
    double frac = pos - index;
    
    if (index + 1 < size) {
        return data_copy[index] + frac * (data_copy[index + 1] - data_copy[index]);
    } else {
        return data_copy[index];
    }
}

/**
 * @brief Calculates the median of an array.
 *
 * @param data Array of double values.
 * @param size Number of elements in the array.
 * @return double The median value.
 */
double calculate_median(double *data, uint16_t size) {
    if (size == 0) return NAN;
    
    // Create a copy to avoid modifying the original data
    double data_copy[MAX_WINDOW_SIZE];
    if (size > MAX_WINDOW_SIZE) size = MAX_WINDOW_SIZE; // Prevent buffer overflow
    memcpy(data_copy, data, sizeof(double) * size);
    
    qsort(data_copy, size, sizeof(double), compare_doubles);
    
    if (size % 2 == 0) {
        return (data_copy[size / 2 - 1] + data_copy[size / 2]) / 2.0;
    } else {
        return data_copy[size / 2];
    }
}

/**
 * @brief Calculates the Median Absolute Deviation (MAD) of an array.
 *
 * @param data Array of double values.
 * @param size Number of elements in the array.
 * @param median The median of the data array.
 * @return double The MAD value. Returns NAN if input is invalid.
 */
double calculate_mad(const double *data, uint16_t size, double median) {
    if (size == 0 || isnan(median)) return NAN;
    
    double deviations[MAX_WINDOW_SIZE];
    uint16_t validCount = 0;
    
    for (uint16_t i = 0; i < size && validCount < MAX_WINDOW_SIZE; ++i) {
        if (isnan(data[i])) continue;
        deviations[validCount++] = fabs(data[i] - median);
    }
    
    if (validCount == 0) return NAN;
    
    // Create a copy to sort
    double deviations_copy[MAX_WINDOW_SIZE];
    memcpy(deviations_copy, deviations, sizeof(double) * validCount);
    
    qsort(deviations_copy, validCount, sizeof(double), compare_doubles);
    
    if (validCount % 2 == 0) {
        return (deviations_copy[validCount / 2 - 1] + deviations_copy[validCount / 2]) / 2.0;
    } else {
        return deviations_copy[validCount / 2];
    }
}

/**
 * @brief Calculates the median and MAD of gradient differences, excluding the top 1% percentile.
 *
 * @param data Array of gradient differences.
 * @param size Number of elements in the data array.
 * @param median Pointer to store the calculated median.
 * @param mad Pointer to store the calculated MAD.
 * @return bool Returns true if calculation was successful, false otherwise.
 */
bool calculateMedianAndMAD(const double *data, uint16_t size, double *median, double *mad) {
    if (!data || size == 0 || !median || !mad) {
        return false;
    }
    
    // Create a copy to sort and manipulate
    double data_copy[MAX_WINDOW_SIZE];
    if (size > MAX_WINDOW_SIZE) size = MAX_WINDOW_SIZE; // Prevent buffer overflow
    memcpy(data_copy, data, sizeof(double) * size);
    
    // Calculate the 99th percentile
    double percentile99 = calculate_percentile(data_copy, size, 99.0);
    
    // Filter out values above the 99th percentile
    double filtered_data[MAX_WINDOW_SIZE];
    uint16_t filtered_size = 0;
    for (uint16_t i = 0; i < size && filtered_size < MAX_WINDOW_SIZE; ++i) {
        if (data[i] <= percentile99) {
            filtered_data[filtered_size++] = data[i];
        }
    }
    
    if (filtered_size == 0) {
        *median = NAN;
        *mad = NAN;
        return false;
    }
    
    // Calculate median of filtered data
    *median = calculate_median(filtered_data, filtered_size);
    
    // Calculate MAD
    *mad = calculate_mad(filtered_data, filtered_size, *median);
    
    return true;
}

/**
 * @file rls_polynomial_regression.h
 * @brief Header for the polynomial regression implementation using RLS and QR decomposition.
 */

#ifndef RLS_POLYNOMIAL_REGRESSION_H
#define RLS_POLYNOMIAL_REGRESSION_H

#include <stdint.h>
#include <stdbool.h>
#include "mqs_def.h"
//#include "mqs/mqs_def.h"

/** Define the size of the sliding window for RLS */
#define RLS_WINDOW 30  /**< Size of the sliding window for RLS */

#define MAX_POLYNOMIAL_DEGREE 3  /**< Maximum degree of polynomial supported (e.g., cubic) */

/**
 * @brief Structure to hold the state of the polynomial regression.
 */
typedef struct {
    uint8_t polynomial_degree;  /**< Degree of the polynomial */
    uint16_t current_num_points;  /**< Current number of data points in the model */
    uint32_t total_data_points_added;  /**< Total number of data points added */
    uint16_t reorthogonalization_counter;  /**< Counter for reorthogonalization */
    double coefficients[MAX_POLYNOMIAL_DEGREE + 1];  /**< Regression coefficients */
    double upper_triangular_R[MAX_POLYNOMIAL_DEGREE + 1][MAX_POLYNOMIAL_DEGREE + 1];  /**< Upper triangular matrix R */
    double Q_transpose_b[MAX_POLYNOMIAL_DEGREE + 1];  /**< Vector Q^T * b */
    int col_permutations[MAX_POLYNOMIAL_DEGREE + 1];  /**< Column permutation indices for pivoting */
    uint16_t oldest_data_index;  /**< Index of the oldest data point in the larger buffer */
    uint16_t max_num_points;  /**< Maximum number of data points (window size) */
} RegressionState;

/**
 * @brief Structure to hold the result of gradient calculations.
 */
typedef struct {
    double gradients[RLS_WINDOW];  /**< Array of calculated gradients */
    uint16_t size;                 /**< Number of gradients calculated */
    double median;                 /**< Median of the gradients */
    double mad;                    /**< Median Absolute Deviation of the gradients */
    uint16_t startIndex;           /**< Start index of the trend in the original data array */
    uint16_t endIndex;             /**< End index of the trend in the original data array */
} GradientCalculationResult;

/**
 * @brief Initializes the RegressionState structure.
 *
 * Allocates memory for the measurement buffer based on the specified window size.
 *
 * @param regression_state Pointer to the RegressionState structure to initialize.
 * @param degree The degree of the polynomial (e.g., 2 for quadratic, 3 for cubic).
 * @param window_size The size of the sliding window for RLS.
 * @return true if initialization is successful, false otherwise.
 */
void initialize_regression_state(RegressionState *regression_state, uint8_t degree, uint16_t max_num_points);

/**
 * @brief Calculates the first-order gradient (slope) of the polynomial function at a specific point.
 *
 * @param regression_state Pointer to the RegressionState structure containing the coefficients.
 * @param x The point at which to calculate the first-order gradient.
 * @return The first-order gradient at the given point x.
 */
double calculate_first_order_gradient(const RegressionState *regression_state, double x);

/**
 * @brief Calculates the second-order gradient (curvature) of the polynomial function at a specific point.
 *
 * @param regression_state Pointer to the RegressionState structure containing the coefficients.
 * @param x The point at which to calculate the second-order gradient.
 * @return The second-order gradient at the given point x.
 */
double calculate_second_order_gradient(const RegressionState *regression_state, double x);

/**
 * @brief Calculates gradients using Recursive Least Squares (RLS) regression.
 *
 * This function iterates through the data points, applies RLS polynomial regression to calculate
 * gradients, and stores them in the `gradCalcResult->gradients` array. It updates the `size`
 * member to reflect the number of gradients calculated.
 *
 * @param values            Array of data points to analyze.
 * @param windowSize        The number of data points to include in the gradient analysis window.
 * @param startIndex        The starting index in the values array.
 * @param degree            The degree of the polynomial used for regression.
 * @param calculate_gradient Function pointer to the gradient calculation function.
 * @param gradCalcResult    Pointer to store the gradient calculation results.
 *
 * @return void
 *
 * @see GradientCalculationResult, calculate_first_order_gradient, calculate_second_order_gradient
 */
void trackGradients(
    const MqsRawDataPoint_t *measurements,
    uint16_t length,
    uint16_t start_index,
    uint8_t degree,
    double (*calculate_gradient)(const RegressionState *, double),
    GradientCalculationResult *result
);



#endif // RLS_POLYNOMIAL_REGRESSION_H

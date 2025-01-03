/**
 * @file rls_polynomial_regression.c
 * @brief Recursive Least Squares (RLS) Polynomial Regression with Hybrid QR Decomposition
 *
 * This file implements a Recursive Least Squares (RLS) algorithm for polynomial regression,
 * optimized for numerical stability and efficiency by combining **Householder reflections**
 * with **Givens rotations** in the QR decomposition process.
 *
 * ## Overview
 *
 * The RLS algorithm updates regression coefficients in real-time as new data points are added
 * or removed. To maintain computational efficiency and numerical accuracy, this implementation
 * employs a hybrid QR decomposition approach:
 *
 * - **Householder Reflections:** Utilized periodically to reorthogonalize the \( \mathbf{Q} \)
 *   matrix, correcting any numerical drift that may accumulate over successive updates.
 *
 * - **Givens Rotations:** Applied incrementally to add or remove data points, allowing
 *   localized adjustments to the upper triangular \( \mathbf{R} \) matrix without
 *   recomputing the entire decomposition.
 *
 * ## Mathematical Foundations
 *
 * ### Householder Reflections
 * - **Nature:** Reflect entire columns to zero out sub-diagonal elements in a single transformation.
 * - **Usage in RLS:** Efficiently reorthogonalizes the \( \mathbf{Q} \) matrix during periodic updates,
 *   ensuring long-term numerical stability of the regression coefficients.
 *
 * ### Givens Rotations
 * - **Nature:** Rotate specific pairs of rows or columns to zero out targeted elements.
 * - **Usage in RLS:** Facilitates incremental updates by efficiently adding or removing data points
 *   through localized transformations, maintaining the upper triangular structure of \( \mathbf{R} \).
 *
 * ## Hybrid QR Decomposition Approach
 *
 * This RLS implementation synergistically combines Householder reflections and Givens rotations
 * to optimize the QR decomposition process, balancing computational efficiency with numerical
 * stability:
 *
 * ### Householder Reflections for Reorthogonalization
 * - **Purpose:** Periodically recompute the entire QR decomposition to correct accumulated numerical
 *   errors and maintain the orthogonality of the \( \mathbf{Q} \) matrix.
 * - **Implementation:** After a set number of incremental updates, Householder reflections are applied
 *   to the current \( \mathbf{Q} \) and \( \mathbf{R} \) matrices to reestablish their orthogonal
 *   and upper triangular properties, respectively.
 *
 * ### Givens Rotations for Incremental Updates
 * - **Purpose:** Efficiently add or remove data points by applying localized rotations, avoiding the
 *   computational overhead of recomputing the entire QR decomposition.
 * - **Implementation:** When a new data point is added, Givens rotations are used to update the
 *   \( \mathbf{R} \) matrix to incorporate the new information. Similarly, when a data point is
 *   removed, Givens rotations adjust \( \mathbf{R} \) to reflect the change, ensuring that
 *   \( \mathbf{Q}\mathbf{R} = \mathbf{A} \) remains valid as \( \mathbf{A} \) evolves.
 *
 * ## Comparative Mathematical Nature
 *
 * | **Aspect**                    | **Householder Reflections**                            | **Givens Rotations**                                   | **Hybrid Approach in RLS**                             |
 * |-------------------------------|--------------------------------------------------------|--------------------------------------------------------|--------------------------------------------------------|
 * | **Transformation Type**      | Reflective transformations on entire columns          | Rotational transformations on specific row/column pairs| Combines global reflections with local rotations       |
 * | **Operation Scope**           | Bulk operations affecting multiple elements simultaneously | Targeted operations affecting specific elements        | Global corrections with localized adjustments          |
 * | **Efficiency**                | Highly efficient for dense, static matrices            | Highly efficient for sparse or dynamic matrices        | Balances efficiency for both dense and dynamic scenarios|
 * | **Numerical Stability**       | High stability for batch processing                    | Maintains stability during incremental updates         | Enhanced stability through periodic reorthogonalization |
 * | **Implementation Complexity** | Simpler for static scenarios                           | More complex due to handling incremental changes        | Increased complexity but optimized for dynamic updates |
 *
 * ## Implications of the Hybrid Approach
 *
 * 1. **Enhanced Efficiency:**
 *    - **Incremental Updates:** Givens rotations enable the QR decomposition to be updated incrementally as new data points are added or old ones are removed, without the need to recompute the entire decomposition.
 *    - **Periodic Reorthogonalization:** Householder reflections periodically reorthogonalize the \( \mathbf{Q} \) matrix, ensuring that numerical errors do not degrade the quality of the decomposition over time.
 *
 * 2. **Improved Numerical Stability:**
 *    - **Combined Strengths:** Givens rotations handle the efficiency of updates, while Householder reflections maintain the numerical accuracy and orthogonality of \( \mathbf{Q} \).
 *    - **Error Mitigation:** Regular reorthogonalization mitigates the accumulation of floating-point errors, preserving the integrity of the regression model.
 *
 * 3. **Flexibility:**
 *    - **Dynamic Data Handling:** The hybrid approach is well-suited for applications where data is continuously streaming in or being removed, such as online learning systems or real-time data analysis.
 *    - **Scalability:** Efficiently handles large datasets by minimizing computational overhead through incremental updates.
 *
 * @note
 * - **Regularization Parameter:** The implementation incorporates a regularization parameter to enhance numerical stability.
 * - **Reorthogonalization Interval:** Determines how frequently Householder reflections are applied to maintain orthogonality.
 * - **Condition Number Monitoring:** Continuously monitors the condition number to detect and address numerical instability.
 *
 * @todo
 * - Implement additional optimization techniques for even greater efficiency.
 * - Extend support for higher-degree polynomials and multi-dimensional regression.
 *
 * @version
 * 1.0
 *
 * @date
 * 2024-11-13
 *
 * @author
 * Tugbars Heptaskin
 */

#include "rls_polynomial_regression.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdarg.h> // For va_list, va_start, va_end

// Uncomment the following line to use global scratch space
#define USE_GLOBAL_SCRATCH_SPACE

/** Regularization parameter for Ridge regression */
#define REGULARIZATION_PARAMETER 1e-4 /**< Regularization parameter lambda */

/** Reorthogonalization interval */
#define REORTHOGONALIZATION_INTERVAL 50 /**< Interval for reorthogonalization */

/** Condition number threshold */
#define CONDITION_NUMBER_THRESHOLD 1e8 /**< Threshold for condition number */

/**
 * @def DEBUG_LEVEL
 * @brief Set the desired debug level (0: No debug, 1: Critical, 2: Detailed, 3: Verbose).
 */
#define DEBUG_LEVEL 0 // Adjust as needed

// Define specific debug print macros
#if DEBUG_LEVEL >= 1
#define DEBUG_PRINT_1(fmt, ...) debug_print(1, fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT_1(fmt, ...) \
    do                          \
    {                           \
    } while (0)
#endif

#if DEBUG_LEVEL >= 2
#define DEBUG_PRINT_2(fmt, ...) debug_print(2, fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT_2(fmt, ...) \
    do                          \
    {                           \
    } while (0)
#endif

#if DEBUG_LEVEL >= 3
#define DEBUG_PRINT_3(fmt, ...) debug_print(3, fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT_3(fmt, ...) \
    do                          \
    {                           \
    } while (0)
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
static void debug_print(int level, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    // Optionally, prepend debug messages with level information
    // printf("[DEBUG LEVEL %d] ", level);
    vprintf(fmt, args);
    va_end(args);
}
#endif

static void recompute_qr_decomposition(
    RegressionState *regression_state,
    const MqsRawDataPoint_t *data);
static inline void solve_for_coefficients(RegressionState *regression_state);
static double compute_condition_number(double R[MAX_POLYNOMIAL_DEGREE + 1][MAX_POLYNOMIAL_DEGREE + 1], int size);
static void updateQR_AddRow(RegressionState *regression_state, const double *new_row, double new_b);
static void updateQR_RemoveRow(RegressionState *regression_state, const double *old_row, double old_b);
void trackFirstOrderGradients(const double *measurements, uint16_t length, uint16_t start_index, uint8_t degree);
void trackSecondOrderGradients(const double *measurements, uint16_t length, uint16_t start_index, uint8_t degree);

#ifdef USE_GLOBAL_SCRATCH_SPACE
// Maximum possible total rows (define this appropriately)
#define MAX_TOTAL_ROWS (RLS_WINDOW + MAX_POLYNOMIAL_DEGREE + 1)

// Declare the arrays as static global variables
static double augmented_matrix_A[MAX_TOTAL_ROWS][MAX_POLYNOMIAL_DEGREE + 1];
static double augmented_vector_b[MAX_TOTAL_ROWS];
#endif

/**
 * @brief Initializes the RegressionState structure.
 *
 * @param regression_state Pointer to the RegressionState structure to initialize.
 * @param degree The degree of the polynomial (e.g., 2 for quadratic, 3 for cubic).
 */
void initialize_regression_state(RegressionState *regression_state, uint8_t degree, uint16_t max_num_points)
{
    DEBUG_PRINT_3("Initializing RegressionState with degree=%u\n", degree);

    regression_state->current_num_points = 0;
    regression_state->max_num_points = max_num_points;
    regression_state->oldest_data_index = 0;
    regression_state->total_data_points_added = 0;
    regression_state->reorthogonalization_counter = 0;
    regression_state->polynomial_degree = degree;

    // Initialize regression coefficients to zero
    memset(regression_state->coefficients, 0, sizeof(regression_state->coefficients));
    DEBUG_PRINT_2("Regression coefficients initialized to zero.\n");

    // Initialize the upper triangular matrix R and Q^T * b vector to zero
    memset(regression_state->upper_triangular_R, 0, sizeof(regression_state->upper_triangular_R));
    memset(regression_state->Q_transpose_b, 0, sizeof(regression_state->Q_transpose_b));
    DEBUG_PRINT_2("Upper triangular matrix R and Q^T * b initialized to zero.\n");

    // Initialize column permutations to default ordering
    for (int i = 0; i <= degree; ++i)
    {
        regression_state->col_permutations[i] = i;
    }
    DEBUG_PRINT_2("Column permutations initialized to default ordering.\n");
}

/**
 * @brief Adds a new data point (phaseAngle) to the RLS regression model and updates its QR decomposition.
 *
 * This function performs an incremental update to the polynomial regression model whenever a new data
 * point becomes available. It uses a QR decomposition-based approach to maintain numerical stability.
 *
 * **Main Steps**:
 *  1. **Prepare the New Data**:
 *     - Extract the measurement (`phaseAngle`) from the provided MqsRawDataPoint_t array.
 *     - Compute the independent variable (`x_value`) based on `total_data_points_added`.
 *     - Build a polynomial basis vector (`new_row`) for the RLS model using `x^0, x^1, ..., x^degree`.
 *
 *  2. **Incremental QR Update**:
 *     - Call `updateQR_AddRow` to incorporate the new row and measurement into the existing QR decomposition.
 *     - If the internal buffer (window) is full, remove the oldest data point using `updateQR_RemoveRow`.
 *
 *  3. **Reorthogonalization**:
 *     - Increment a counter each time we add a point. Once a threshold (`REORTHOGONALIZATION_INTERVAL`)
 *       is reached, recompute the entire QR decomposition from scratch for higher accuracy.
 *
 *  4. **Regression Coefficients**:
 *     - Solve for the new regression coefficients via `solve_for_coefficients`.
 *     - Monitor the condition number to detect potential numerical instability. If it exceeds a threshold,
 *       recompute the QR decomposition again and re-solve the coefficients.
 *
 * **Usage**:
 * This function is typically called whenever a new measurement arrives (e.g., in a real-time system).
 * If your buffer is not yet full, it simply increases `current_num_points` until it reaches `max_num_points`.
 * Once full, it removes the oldest data point to maintain a sliding window approach.
 *
 * @param[in,out] regression_state
 *    Pointer to the `RegressionState` structure that stores the current QR decomposition, polynomial degree,
 *    buffer sizes, and other relevant state. This state is modified in place.
 *
 * @param[in] data
 *    Pointer to an array of `MqsRawDataPoint_t`, from which we extract the `.phaseAngle` value for the new
 *    measurement. The array must be at least as large as the largest `data_index` you will use.
 *
 * @param[in] data_index
 *    The index in the `data` array from which to read the new measurement's `.phaseAngle`.
 *    This identifies which `phaseAngle` to insert into the regression model.
 */
void add_data_point_to_regression(
    RegressionState *regression_state,
    const MqsRawDataPoint_t *data,
    uint16_t data_index)
{
    // (A) Determine how many polynomial coefficients we need (polynomial_degree + 1).
    int num_coefficients = regression_state->polynomial_degree + 1;

    // (B) Extract the new measurement from data[data_index].phaseAngle.
    double measurement = (double)data[data_index].phaseAngle;
    DEBUG_PRINT_3("Adding new data point: measurement=%.6f at data_index=%u\n", measurement, data_index);

    // (C) Compute the x_value for the new measurement (often just the index in time).
    double x_value = (double)(regression_state->total_data_points_added);
    DEBUG_PRINT_2("Current x_value for new data point: %.2f\n", x_value);

    // (D) Build the polynomial basis vector (e.g., [1, x, x^2, ... x^degree]).
    double new_row[MAX_POLYNOMIAL_DEGREE + 1];
    double x_power = 1.0;
    for (int i = 0; i < num_coefficients; ++i)
    {
        new_row[i] = x_power;
        x_power *= x_value;
    }
    DEBUG_PRINT_3("New polynomial basis vector computed.\n");

    // (E) Add the new row and measurement to the existing QR decomposition incrementally.
    updateQR_AddRow(regression_state, new_row, measurement);
    DEBUG_PRINT_3("QR decomposition updated with new data point.\n");

    // (F) Check if we've exceeded the maximum number of points in our sliding window.
    if (regression_state->current_num_points >= regression_state->max_num_points)
    {
        // (F.1) We must remove the oldest data point from the model to maintain the window size.
        uint16_t oldest_data_index = regression_state->oldest_data_index;
        double old_measurement = (double)data[oldest_data_index].phaseAngle;
        DEBUG_PRINT_3("Buffer full. Removing oldest data point at index %u with measurement=%.6f.\n",
                      oldest_data_index, old_measurement);

        // (F.2) Compute the polynomial basis for the oldest data point's x_value.
        double old_x_value = (double)(regression_state->total_data_points_added - regression_state->max_num_points);
        double old_row[MAX_POLYNOMIAL_DEGREE + 1];
        x_power = 1.0;
        for (int i = 0; i < num_coefficients; ++i)
        {
            old_row[i] = x_power;
            x_power *= old_x_value;
        }
        DEBUG_PRINT_3("Old polynomial basis vector computed for removal.\n");

        // (F.3) Perform a downdate on the QR decomposition to remove the effect of that old data point.
        updateQR_RemoveRow(regression_state, old_row, old_measurement);
        DEBUG_PRINT_3("QR decomposition updated by removing oldest data point.\n");

        // (F.4) Advance the oldest data index, so we know which point to remove next time.
        regression_state->oldest_data_index = (uint16_t)(oldest_data_index + 1);
        DEBUG_PRINT_3("Oldest data index updated to %u.\n", regression_state->oldest_data_index);
    }
    else
    {
        // (G) If we're not yet at the max number of points, just increment current_num_points.
        regression_state->current_num_points++;
        DEBUG_PRINT_3("Buffer not full yet. Current number of points: %u.\n",
                      regression_state->current_num_points);
    }

    // (H) We've successfully added a new data point, so increment total_data_points_added.
    regression_state->total_data_points_added++;
    DEBUG_PRINT_3("Total data points added incremented to %u.\n",
                  regression_state->total_data_points_added);

    // (I) Consider whether we need to reorthogonalize from scratch. We do so every REORTHOGONALIZATION_INTERVAL points.
    regression_state->reorthogonalization_counter++;
    DEBUG_PRINT_3("Reorthogonalization counter incremented to %u.\n",
                  regression_state->reorthogonalization_counter);

    if (regression_state->reorthogonalization_counter >= REORTHOGONALIZATION_INTERVAL)
    {
        DEBUG_PRINT_3("Reorthogonalization interval reached. Recomputing QR decomposition from scratch.\n");
        // Rebuild the entire QR decomposition for better numerical stability.
        recompute_qr_decomposition(regression_state, data);
        regression_state->reorthogonalization_counter = 0;
        DEBUG_PRINT_3("QR decomposition recomputed and counter reset.\n");
    }

    // (J) Solve for the new polynomial regression coefficients now that the model is updated.
    solve_for_coefficients(regression_state);
    DEBUG_PRINT_2("Regression coefficients solved and updated.\n");

    // (K) Compute the condition number to monitor numerical stability.
    double condition_number = compute_condition_number(
        regression_state->upper_triangular_R,
        num_coefficients);
    DEBUG_PRINT_2("Computed condition number: %.6e\n", condition_number);

    // (L) If the condition number is too large, we rebuild the entire decomposition and re-solve.
    if (condition_number > CONDITION_NUMBER_THRESHOLD)
    {
        DEBUG_PRINT_1("Condition number %.6e exceeds threshold %.6e. Recomputing QR.\n",
                      condition_number, CONDITION_NUMBER_THRESHOLD);

        // Redo the entire QR from scratch
        recompute_qr_decomposition(regression_state, data);
        solve_for_coefficients(regression_state);
        DEBUG_PRINT_2("Regression coefficients re-solved after recomputing QR decomposition.\n");
    }
}

/**
 * @brief Updates the QR decomposition by adding a new row using Givens rotations.
 *
 * This function incrementally updates the existing QR decomposition when a new data point is added
 * to the regression model. It uses Givens rotations to maintain the upper triangular
 * structure of the R matrix without recomputing the entire decomposition.
 *
 * **Changes Introduced:**
 * - **Use of Givens Rotations:**
 *   - Previously, the entire QR decomposition was recomputed when new data was added.
 *   - Now, Givens rotations are employed to update the QR factors incrementally.
 *   - This change improves efficiency by reducing computational load and enables real-time updates.
 * - **Efficiency and Numerical Stability:**
 *   - Givens rotations are numerically stable for incremental updates, especially when dealing
 *     with streaming data or large datasets.
 *
 * **Why Givens Rotations Were Introduced:**
 * - To allow efficient incremental updates to the QR decomposition without full recomputation.
 * - To enhance numerical stability during the addition of new data points.
 *
 * @param regression_state Pointer to the RegressionState structure.
 * @param new_row The new row to add (polynomial basis vector).
 * @param new_b The new measurement value.
 */
static void updateQR_AddRow(RegressionState *regression_state, const double *new_row, double new_b)
{
    int num_coefficients = regression_state->polynomial_degree + 1;
    DEBUG_PRINT_3("Updating QR decomposition by adding a new row.\n");

    // Copy the new row into a temporary array
    double r_row[MAX_POLYNOMIAL_DEGREE + 1];
    memcpy(r_row, new_row, sizeof(double) * num_coefficients);
    DEBUG_PRINT_2("New row copied for QR update.\n");

    double Q_tb_new = new_b;
    DEBUG_PRINT_2("Initial Q^T * b for new row: %.6f\n", Q_tb_new);

    // Apply Givens rotations to zero out sub-diagonal elements
    for (int i = 0; i < num_coefficients; ++i)
    {
        double a = regression_state->upper_triangular_R[i][i];
        double b = r_row[i];

        if (fabs(b) > 1e-10)
        { // Avoid division by zero
            // Compute Givens rotation parameters
            double r = hypot(a, b);
            double c = a / r;
            double s = b / r;

            DEBUG_PRINT_2("Applying Givens rotation at column %d: c=%.6f, s=%.6f\n", i, c, s);

            // Update R
            regression_state->upper_triangular_R[i][i] = r;
            for (int j = i + 1; j < num_coefficients; ++j)
            {
                double temp = regression_state->upper_triangular_R[i][j];
                regression_state->upper_triangular_R[i][j] = c * temp + s * r_row[j];
                r_row[j] = -s * temp + c * r_row[j];
                DEBUG_PRINT_3("Updated R[%d][%d]=%.6f and r_row[%d]=%.6f\n", i, j, regression_state->upper_triangular_R[i][j], j, r_row[j]);
            }

            // Update Q^T * b
            double temp_b = regression_state->Q_transpose_b[i];
            regression_state->Q_transpose_b[i] = c * temp_b + s * Q_tb_new;
            Q_tb_new = -s * temp_b + c * Q_tb_new;
            DEBUG_PRINT_3("Updated Q_transpose_b[%d]=%.6f and Q_tb_new=%.6f\n", i, regression_state->Q_transpose_b[i], Q_tb_new);
        }
        else
        {
            // No rotation needed; R[i][i] remains the same
            DEBUG_PRINT_3("No rotation needed for column %d.\n", i);
            // Update Q^T * b with the new data
            regression_state->Q_transpose_b[i] += regression_state->upper_triangular_R[i][i] * 0.0;
        }
    }

    DEBUG_PRINT_2("QR decomposition updated with new row.\n");
}

/**
 * @brief Updates the QR decomposition by removing an old row using Givens rotations.
 *
 * This function incrementally downdates the existing QR decomposition when an old data point
 * is removed from the regression model. It uses Givens rotations to adjust the R matrix,
 * maintaining its upper triangular structure.
 *
 * **Changes Introduced:**
 * - **Use of Givens Rotations for Downdating:**
 *   - Previously, the QR decomposition was recomputed entirely when old data was removed.
 *   - Now, Givens rotations are used to efficiently remove the influence of the oldest data point.
 * - **Efficiency Improvement:**
 *   - This change reduces computational overhead, making the algorithm suitable for real-time applications.
 * - **Numerical Stability:**
 *   - Givens rotations provide a stable method for downdating, minimizing numerical errors.
 *
 * **Why Givens Rotations Were Introduced:**
 * - To enable efficient and stable downdating of the QR decomposition when data points are removed.
 * - To maintain the numerical integrity of the regression model over time.
 *
 * @param regression_state Pointer to the RegressionState structure.
 * @param old_row The old row to remove (polynomial basis vector).
 * @param old_b The old measurement value.
 */
static void updateQR_RemoveRow(RegressionState *regression_state, const double *old_row, double old_b)
{
    int num_coefficients = regression_state->polynomial_degree + 1;
    DEBUG_PRINT_3("Updating QR decomposition by removing an old row.\n");

    // Subtract the influence of the old data point
    double r_row[MAX_POLYNOMIAL_DEGREE + 1];
    memcpy(r_row, old_row, sizeof(double) * num_coefficients);
    DEBUG_PRINT_2("Old row copied for QR downdate.\n");

    double Q_tb_old = old_b;
    DEBUG_PRINT_2("Initial Q^T * b for old row: %.6f\n", Q_tb_old);

    // Apply reverse Givens rotations to remove the old data point
    for (int i = 0; i < num_coefficients; ++i)
    {
        double a = regression_state->upper_triangular_R[i][i];
        double b = r_row[i];

        if (fabs(b) > 1e-10)
        { // Avoid division by zero
            // Compute Givens rotation parameters
            double r = hypot(a, b);
            double c = a / r;
            double s = b / r;

            DEBUG_PRINT_2("Applying Givens rotation at column %d: c=%.6f, s=%.6f\n", i, c, s);

            // Update R
            regression_state->upper_triangular_R[i][i] = r;
            for (int j = i + 1; j < num_coefficients; ++j)
            {
                double temp = regression_state->upper_triangular_R[i][j];
                regression_state->upper_triangular_R[i][j] = c * temp - s * r_row[j];
                r_row[j] = s * temp + c * r_row[j];
                DEBUG_PRINT_3("Updated R[%d][%d]=%.6f and r_row[%d]=%.6f\n", i, j, regression_state->upper_triangular_R[i][j], j, r_row[j]);
            }

            // Update Q^T * b
            double temp_b = regression_state->Q_transpose_b[i];
            regression_state->Q_transpose_b[i] = c * temp_b - s * Q_tb_old;
            Q_tb_old = s * temp_b + c * Q_tb_old;
            DEBUG_PRINT_3("Updated Q_transpose_b[%d]=%.6f and Q_tb_old=%.6f\n", i, regression_state->Q_transpose_b[i], Q_tb_old);
        }
        else
        {
            // No rotation needed
            DEBUG_PRINT_3("No rotation needed for column %d.\n", i);
            // Update Q^T * b
            // regression_state->Q_transpose_b[i] -= regression_state->upper_triangular_R[i][i] * 0.0;
        }
    }

    DEBUG_PRINT_2("QR decomposition updated by removing old row.\n");
}

/**
 * @brief Recomputes the QR decomposition using optimized Householder reflections with column pivoting.
 *
 * This function performs a full recomputation of the QR decomposition using Householder reflections
 * with column pivoting for enhanced numerical stability. Column pivoting rearranges the columns
 * based on their norms to position the largest possible pivot elements on the diagonal, reducing
 * potential numerical issues.
 *
 * **Changes Introduced**:
 * - **Access to External Phase-Angle Data**:
 *   - The function now accepts a `MqsRawDataPoint_t *data` array instead of a `double *`.
 *   - The measurement (dependent variable) is extracted from `data[i].phaseAngle`.
 * - **Column Pivoting**:
 *   - Implements column pivoting by tracking and swapping columns based on their norms.
 *   - Updates the column indices array in `RegressionState` to keep track of column permutations.
 * - **Regularization and Minimizing Stack Usage**:
 *   - Adds diagonal regularization rows to mitigate overfitting.
 *   - Optionally uses global scratch space (`USE_GLOBAL_SCRATCH_SPACE`) to avoid large local arrays.
 *
 * **Why These Changes Were Introduced**:
 * - To eliminate redundant storage of measurements within `RegressionState`.
 * - To handle real-time updates from an external buffer of MqsRawDataPoint_t (i.e., `.phaseAngle`).
 * - To enhance numerical stability during QR decomposition (column pivoting).
 * - To minimize stack usage by optionally using global arrays.
 *
 * **Algorithmic Steps**:
 *  1. **Setup / Data Extraction**:
 *     - Reads `num_data_points` from `regression_state->current_num_points`.
 *     - Each data point’s index is computed as `oldest_data_index + i`.
 *     - The independent variable `x` is computed as `(total_data_points_added - num_data_points + i)`.
 *     - The dependent variable `measurement` is read from `data[data_index].phaseAngle`.
 *  2. **Construct Augmented Matrix**:
 *     - Builds a polynomial design matrix \( A \) of size \((\text{num_data_points} + \text{num_coefficients}) \times \text{num_coefficients}\).
 *     - Appends `num_coefficients` regularization rows to reduce overfitting.
 *  3. **Initialize QR**:
 *     - Resets `upper_triangular_R` and `Q_transpose_b` in `RegressionState`.
 *     - Prepares column norms and an array of column indices for pivoting.
 *  4. **Householder Reflections (with Pivoting)**:
 *     - For each column \( k \):
 *       - Find the column of max norm among [k..end].
 *       - Swap columns (and norms, col_indices) if necessary.
 *       - Form and apply the Householder transformation to zero out below-diagonal elements.
 *       - Update the augmented vector \( b \) (i.e., `augmented_vector_b`) accordingly.
 *  5. **Extract R, Qᵀb, and Column Permutations**:
 *     - Copies the upper triangular portion of \( A \) into `upper_triangular_R`.
 *     - Copies the transformed vector into `Q_transpose_b`.
 *     - Stores the pivoting permutation array in `regression_state->col_permutations`.
 *
 * **Result**:
 * - An updated QR factorization is available in `regression_state->upper_triangular_R` and
 *   `regression_state->Q_transpose_b`, along with column permutations for correct coefficient mapping.
 *
 * @param[in,out] regression_state
 *    Pointer to the `RegressionState` structure holding the current polynomial degree,
 *    buffer indices, and storage for R and Qᵀb. This function updates those fields.
 *
 * @param[in] data
 *    External buffer of `MqsRawDataPoint_t`. The measurement (dependent variable) is taken
 *    from `data[data_index].phaseAngle`. The order of these data points corresponds to
 *    the window from `oldest_data_index` up to `current_num_points`.
 */
static void recompute_qr_decomposition(
    RegressionState *regression_state,
    const MqsRawDataPoint_t *data)
{
    int num_data_points = regression_state->current_num_points;
    int num_coefficients = regression_state->polynomial_degree + 1;
    DEBUG_PRINT_3("Recomputing QR decomposition with column pivoting for RegressionState.\n");

    // Regularization parameter to prevent overfitting
    double regularization_lambda = REGULARIZATION_PARAMETER;
    DEBUG_PRINT_2("Regularization parameter lambda=%.6f\n", regularization_lambda);

    // Total number of rows after adding regularization terms
    int total_rows = num_data_points + num_coefficients;

    // Maximum possible total rows
    // int max_total_rows = regression_state->max_num_points + num_coefficients;

#ifdef USE_GLOBAL_SCRATCH_SPACE
    // Use global scratch arrays: augmented_matrix_A, augmented_vector_b
    // Validate that total_rows does not exceed our pre-allocated maximum
    if (total_rows > MAX_TOTAL_ROWS)
    {
        DEBUG_PRINT_1("Error: total_rows (%d) exceeds MAX_TOTAL_ROWS (%d).\n", total_rows, MAX_TOTAL_ROWS);
        return;
    }
#else
    // Local stack-based arrays for the augmented matrix and vector
    double augmented_matrix_A[total_rows][MAX_POLYNOMIAL_DEGREE + 1];
    double augmented_vector_b[total_rows];
#endif

    // Initialize the augmented arrays to zero
    memset(augmented_matrix_A, 0, sizeof(double) * total_rows * (MAX_POLYNOMIAL_DEGREE + 1));
    memset(augmented_vector_b, 0, sizeof(double) * total_rows);
    DEBUG_PRINT_2("Augmented matrix A and vector b initialized to zero.\n");

    // Initialize column indices for pivoting
    int col_indices[MAX_POLYNOMIAL_DEGREE + 1];
    for (int i = 0; i < num_coefficients; ++i)
    {
        col_indices[i] = i;
    }

    // Populate the augmented matrix with the data from [oldest_data_index ... oldest_data_index + num_data_points - 1]
    for (int i = 0; i < num_data_points; ++i)
    {
        uint16_t data_index = (uint16_t)(regression_state->oldest_data_index + i);

        // The x_value is computed based on how many total points have been added so far,
        // minus how many points we have, plus i
        double x = (double)(regression_state->total_data_points_added - num_data_points + i);

        // Extract the measurement from data[data_index].phaseAngle
        double measurement = (double)data[data_index].phaseAngle;
        DEBUG_PRINT_3("Processing data point i=%d: x=%.2f, measurement=%.6f (data_index=%u)\n",
                      i, x, measurement, data_index);

        // Fill row i of A with polynomial terms 1, x, x^2, ...
        double x_power = 1.0;
        for (int j = 0; j < num_coefficients; ++j)
        {
            augmented_matrix_A[i][j] = x_power;
            x_power *= x;
        }

        // Put the measurement in augmented_vector_b
        augmented_vector_b[i] = measurement;
    }
    DEBUG_PRINT_2("Augmented matrix populated with original data.\n");

    // Add regularization rows at the bottom of the matrix
    for (int i = 0; i < num_coefficients; ++i)
    {
        int row = num_data_points + i;
        // sqrt(lambda) * I
        augmented_matrix_A[row][i] = sqrt(regularization_lambda);
        augmented_vector_b[row] = 0.0;
        DEBUG_PRINT_3("Added regularization row %d: A[%d][%d]=%.6f, b[%d]=%.6f\n",
                      i, row, i, augmented_matrix_A[row][i], row, augmented_vector_b[row]);
    }
    DEBUG_PRINT_2("Regularization rows added.\n");

    // Reset R and Qᵀb in the regression state
    memset(regression_state->upper_triangular_R, 0, sizeof(regression_state->upper_triangular_R));
    memset(regression_state->Q_transpose_b, 0, sizeof(regression_state->Q_transpose_b));
    DEBUG_PRINT_2("Reset upper_triangular_R and Q_transpose_b.\n");

    // Compute column norms for pivoting
    double col_norms[MAX_POLYNOMIAL_DEGREE + 1];
    for (int j = 0; j < num_coefficients; ++j)
    {
        double sum = 0.0;
        for (int i = 0; i < total_rows; ++i)
        {
            sum += augmented_matrix_A[i][j] * augmented_matrix_A[i][j];
        }
        col_norms[j] = sqrt(sum);
    }

    // Householder QR with Column Pivoting
    for (int k = 0; k < num_coefficients; ++k)
    {
        // 1) Find the column with the maximum norm among [k..(num_coefficients-1)]
        int max_col = k;
        double max_norm = col_norms[k];
        for (int j = k + 1; j < num_coefficients; ++j)
        {
            if (col_norms[j] > max_norm)
            {
                max_norm = col_norms[j];
                max_col = j;
            }
        }

        // 2) Swap columns in augmented_matrix_A, col_indices, and col_norms if needed
        if (max_col != k)
        {
            for (int i = 0; i < total_rows; ++i)
            {
                double temp = augmented_matrix_A[i][k];
                augmented_matrix_A[i][k] = augmented_matrix_A[i][max_col];
                augmented_matrix_A[i][max_col] = temp;
            }
            // col indices
            int temp_idx = col_indices[k];
            col_indices[k] = col_indices[max_col];
            col_indices[max_col] = temp_idx;
            // norms
            double temp_norm = col_norms[k];
            col_norms[k] = col_norms[max_col];
            col_norms[max_col] = temp_norm;

            DEBUG_PRINT_2("Swapped columns %d and %d for pivoting.\n", k, max_col);
        }

        // 3) Householder reflect on column k to zero out below the diagonal
        double sigma = 0.0;
        for (int i = k; i < total_rows; ++i)
        {
            double val = augmented_matrix_A[i][k];
            sigma += val * val;
        }
        sigma = sqrt(sigma);
        DEBUG_PRINT_3("Householder transformation at column %d: sigma=%.6f\n", k, sigma);

        if (sigma == 0.0)
        {
            DEBUG_PRINT_3("Column %d has zero norm. Skipping.\n", k);
            continue;
        }

        // Form the Householder vector
        double vk = augmented_matrix_A[k][k] + (augmented_matrix_A[k][k] >= 0.0 ? sigma : -sigma);

        if (fabs(vk) < 1e-12)
        {
            DEBUG_PRINT_3("vk too small at col %d. Skipping transformation.\n", k);
            continue;
        }

        double beta = 1.0 / (sigma * vk);

        // Apply to the column
        for (int i = k + 1; i < total_rows; ++i)
        {
            augmented_matrix_A[i][k] /= vk;
        }
        augmented_matrix_A[k][k] = -sigma;

        // Apply transformation to subsequent columns
        for (int j = k + 1; j < num_coefficients; ++j)
        {
            double s = 0.0;
            for (int i = k; i < total_rows; ++i)
            {
                s += augmented_matrix_A[i][k] * augmented_matrix_A[i][j];
            }
            s *= beta;
            for (int i = k; i < total_rows; ++i)
            {
                augmented_matrix_A[i][j] -= augmented_matrix_A[i][k] * s;
            }
        }

        // Apply transformation to vector b
        double s = 0.0;
        for (int i = k; i < total_rows; ++i)
        {
            s += augmented_matrix_A[i][k] * augmented_vector_b[i];
        }
        s *= beta;
        for (int i = k; i < total_rows; ++i)
        {
            augmented_vector_b[i] -= augmented_matrix_A[i][k] * s;
        }

        // Zero out below-diagonal elements in column k
        for (int i = k + 1; i < total_rows; ++i)
        {
            augmented_matrix_A[i][k] = 0.0;
        }

        // Update norms for remaining columns
        for (int j = k + 1; j < num_coefficients; ++j)
        {
            double sum = 0.0;
            for (int i = k + 1; i < total_rows; ++i)
            {
                sum += augmented_matrix_A[i][j] * augmented_matrix_A[i][j];
            }
            col_norms[j] = sqrt(sum);
        }
    }

    // 4) Extract R and Qᵀb from the transformed matrix
    for (int i = 0; i < num_coefficients; ++i)
    {
        for (int j = i; j < num_coefficients; ++j)
        {
            regression_state->upper_triangular_R[i][j] = augmented_matrix_A[i][j];
        }
        regression_state->Q_transpose_b[i] = augmented_vector_b[i];
    }

    // 5) Store the column permutation indices in the regression state
    memcpy(regression_state->col_permutations, col_indices,
           sizeof(int) * num_coefficients);

    DEBUG_PRINT_3("QR decomposition recomputed with column pivoting (using MqsRawDataPoint_t).\n");
}

/**
 * @brief Solves the upper triangular system to find the regression coefficients with column permutations.
 *
 * This function performs back substitution on the upper triangular matrix R to solve for
 * the regression coefficients, taking into account any column permutations due to pivoting.
 * It ensures that the latest coefficients are used for gradient calculations and predictions.
 *
 * **Changes Introduced:**
 * - **Handling Column Permutations:**
 *   - Adjusts the back substitution process to account for the column permutations.
 *   - Rearranges the coefficients to match the original variable ordering.
 *
 * **Why These Changes Were Introduced:**
 * - To correctly solve for coefficients when column pivoting is used.
 * - To ensure that the coefficients correspond to the correct polynomial terms.
 *
 * @param regression_state Pointer to the RegressionState structure.
 */
static inline void solve_for_coefficients(RegressionState *regression_state)
{
    int num_coefficients = regression_state->polynomial_degree + 1;
    DEBUG_PRINT_3("Solving for regression coefficients using back substitution with column permutations.\n");

    double temp_coefficients[MAX_POLYNOMIAL_DEGREE + 1];

    // Perform back substitution to solve R * x = Q^T * b
    for (int i = num_coefficients - 1; i >= 0; --i)
    {
        double sum = regression_state->Q_transpose_b[i];
        for (int j = i + 1; j < num_coefficients; ++j)
        {
            sum -= regression_state->upper_triangular_R[i][j] * temp_coefficients[j];
        }

        if (fabs(regression_state->upper_triangular_R[i][i]) < 1e-12)
        {
            DEBUG_PRINT_1("Warning: Small diagonal element R[%d][%d]=%.6e during back substitution.\n", i, i, regression_state->upper_triangular_R[i][i]);
            temp_coefficients[i] = 0.0; // Assign zero to avoid division by zero
        }
        else
        {
            temp_coefficients[i] = sum / regression_state->upper_triangular_R[i][i];
            DEBUG_PRINT_3("Temporary Coefficient[%d] solved: %.6f\n", i, temp_coefficients[i]);
        }
    }

    // Rearrange the solution according to the original column order
    for (int i = 0; i < num_coefficients; ++i)
    {
        int permuted_index = regression_state->col_permutations[i];
        regression_state->coefficients[permuted_index] = temp_coefficients[i];
        DEBUG_PRINT_3("Coefficient[%d] set to %.6f (permuted from position %d)\n", permuted_index, temp_coefficients[i], i);
    }

    DEBUG_PRINT_2("Regression coefficients updated with column permutations.\n");
}

#ifdef USE_IMPROVED_CONDITION_NUMBER
/**
 * @brief Computes the 1-norm of the upper triangular matrix R.
 *
 * This function calculates the maximum absolute column sum of the upper triangular matrix R,
 * which is useful for estimating the condition number of R.
 *
 * **Why This Function Was Introduced:**
 * - To provide a more accurate estimation of the matrix norm, which is essential for computing the condition number.
 * - Utilizing the 1-norm leverages the structure of the upper triangular matrix for efficient computation.
 *
 * @param R The upper triangular matrix R.
 * @param size The size of the matrix (number of coefficients).
 * @return The 1-norm of matrix R.
 */
static double compute_R_norm_1(const double R[MAX_POLYNOMIAL_DEGREE + 1][MAX_POLYNOMIAL_DEGREE + 1], int size)
{
    DEBUG_PRINT_3("Computing the 1-norm of matrix R.\n");
    double norm = 0.0;
    for (int j = 0; j < size; ++j)
    {
        double sum = 0.0;
        for (int i = 0; i <= j; ++i)
        { // Only upper triangular part
            sum += fabs(R[i][j]);
            DEBUG_PRINT_3("Adding |R[%d][%d]|=%.6f to column sum.\n", i, j, fabs(R[i][j]));
        }
        DEBUG_PRINT_3("Column %d sum: %.6f\n", j, sum);
        if (sum > norm)
        {
            norm = sum;
            DEBUG_PRINT_3("Updated norm to %.6f\n", norm);
        }
    }
    DEBUG_PRINT_2("Computed 1-norm of R: %.6f\n", norm);
    return norm;
}

/**
 * @brief Estimates the 1-norm of the inverse of the upper triangular matrix R.
 *
 * This function estimates the maximum absolute column sum of the inverse of R without explicitly computing R^{-1}.
 * It uses an iterative method suitable for upper triangular matrices.
 *
 * **Why This Function Was Introduced:**
 * - To enable accurate estimation of the condition number by computing \|R^{-1}\|_1 efficiently.
 * - Avoids explicit inversion of R, which can be computationally expensive and numerically unstable.
 *
 * @param R The upper triangular matrix R.
 * @param size The size of the matrix (number of coefficients).
 * @return The estimated 1-norm of R^{-1}.
 */
static double estimate_R_inverse_norm_1(const double R[MAX_POLYNOMIAL_DEGREE + 1][MAX_POLYNOMIAL_DEGREE + 1], int size)
{
    DEBUG_PRINT_3("Estimating the 1-norm of the inverse of matrix R.\n");
    double x[MAX_POLYNOMIAL_DEGREE + 1];
    double z[MAX_POLYNOMIAL_DEGREE + 1];
    double est = 0.0;

    // Initialize vectors
    for (int i = 0; i < size; ++i)
    {
        x[i] = 1.0 / (double)size;
        DEBUG_PRINT_3("Initialized x[%d]=%.6f\n", i, x[i]);
    }

    int iter = 0;
    double prev_est = 0.0;
    const int max_iter = 5;
    const double tol = 1e-6;

    do
    {
        // Solve R * z = x using back substitution
        for (int i = size - 1; i >= 0; --i)
        {
            double sum = x[i];
            for (int j = i + 1; j < size; ++j)
            {
                sum -= R[i][j] * z[j];
                DEBUG_PRINT_3("Subtracting R[%d][%d]*z[%d]=%.6f*%.6f from sum.\n", i, j, j, R[i][j], z[j]);
            }
            z[i] = sum / R[i][i];
            DEBUG_PRINT_3("Computed z[%d]=%.6f\n", i, z[i]);
        }

        // Compute the 1-norm of z
        est = 0.0;
        for (int i = 0; i < size; ++i)
        {
            est += fabs(z[i]);
            DEBUG_PRINT_3("Adding |z[%d]|=%.6f to est.\n", i, fabs(z[i]));
        }
        DEBUG_PRINT_2("Current estimate of norm_R_inv: %.6f\n", est);

        // Normalize z
        double z_norm = est;
        for (int i = 0; i < size; ++i)
        {
            z[i] /= z_norm;
            DEBUG_PRINT_3("Normalized z[%d]=%.6f\n", i, z[i]);
        }

        // Check for convergence
        if (fabs(est - prev_est) < tol * est)
        {
            DEBUG_PRINT_2("Convergence achieved after %d iterations.\n", iter);
            break;
        }

        // Prepare for next iteration
        memcpy(x, z, sizeof(double) * size);
        prev_est = est;
        iter++;
        DEBUG_PRINT_3("Iteration %d completed. est=%.6f, prev_est=%.6f\n", iter, est, prev_est);
    } while (iter < max_iter);

    DEBUG_PRINT_2("Estimated 1-norm of R^{-1}: %.6f after %d iterations.\n", est, iter);
    return est;
}

/**
 * @brief Computes the condition number of the upper triangular matrix R using the 1-norm.
 *
 * This function computes the condition number κ(R) = \|R\|₁ * \|R⁻¹\|₁, providing
 * a more accurate estimation compared to using only the diagonal elements.
 *
 * **Changes Introduced:**
 * - Uses matrix norms to estimate the condition number more accurately.
 * - Introduces efficient estimation of \|R⁻¹\|₁ without explicit inversion.
 *
 * **Why These Changes Were Introduced:**
 * - To improve the reliability of the condition number estimation.
 * - Enhances the ability to detect numerical instability in the regression model.
 *
 * @param R The upper triangular matrix R.
 * @param size The size of the matrix (number of coefficients).
 * @return The estimated condition number.
 */
static double compute_condition_number(const double R[MAX_POLYNOMIAL_DEGREE + 1][MAX_POLYNOMIAL_DEGREE + 1], int size)
{
    DEBUG_PRINT_3("Computing condition number of matrix R using 1-norm.\n");
    double norm_R = compute_R_norm_1(R, size);
    double norm_R_inv = estimate_R_inverse_norm_1(R, size);
    double condition_number = norm_R * norm_R_inv;
    DEBUG_PRINT_2("Computed condition number: %.6e (norm_R=%.6e, norm_R_inv=%.6e)\n", condition_number, norm_R, norm_R_inv);
    return condition_number;
}
#else  // Original method
/**
 * @brief Computes the condition number of the upper triangular matrix R using the 2-norm estimation.
 *
 * @param R The upper triangular matrix R.
 * @param size The size of the matrix (number of coefficients).
 * @return The estimated condition number.
 */
static double compute_condition_number(double R[MAX_POLYNOMIAL_DEGREE + 1][MAX_POLYNOMIAL_DEGREE + 1], int size)
{
    DEBUG_PRINT_3("Computing condition number of matrix R.\n");

    // Estimate the condition number using the ratio of norms
    // Since R is upper triangular, we can estimate the condition number by the ratio of largest to smallest singular values
    // Use the absolute values of diagonal elements as approximations
    double max_sv = 0.0;
    double min_sv = DBL_MAX;

    for (int i = 0; i < size; ++i)
    {
        double diag = fabs(R[i][i]);
        if (diag > max_sv)
        {
            max_sv = diag;
        }
        if (diag < min_sv && diag > 1e-12)
        {
            min_sv = diag;
        }
    }

    DEBUG_PRINT_3("Maximum singular value estimate: %.6e, Minimum singular value estimate: %.6e\n", max_sv, min_sv);

    if (min_sv < 1e-12)
    {
        // Matrix is singular or nearly singular
        DEBUG_PRINT_1("Condition number is infinite (matrix is singular or nearly singular).\n");
        return DBL_MAX;
    }

    double condition_number = max_sv / min_sv;
    DEBUG_PRINT_3("Computed condition number: %.6e\n", condition_number);

    return condition_number;
}
#endif // USE_IMPROVED_CONDITION_NUMBER

/**
 * @brief Calculates the first-order gradient (slope) of the polynomial function at a specific point.
 *
 * @param regression_state Pointer to the RegressionState structure containing the coefficients.
 * @param x The point at which to calculate the first-order gradient.
 * @return The first-order gradient at the given point x.
 */
double calculate_first_order_gradient(const RegressionState *regression_state, double x)
{
    DEBUG_PRINT_3("Calculating first-order gradient at x=%.6f\n", x);

    double derivative = 0.0;
    double x_power = 1.0;
    int degree = regression_state->polynomial_degree;

    for (int i = 1; i <= degree; ++i)
    {
        x_power *= x;
        derivative += i * regression_state->coefficients[i] * x_power / x;
        DEBUG_PRINT_3("Term %d: i=%d, coefficient=%.6f, x_power=%.6f, derivative=%.6f\n", i, i, regression_state->coefficients[i], x_power, derivative);
    }

    DEBUG_PRINT_3("First-order gradient calculated: %.6f\n", derivative);
    return derivative;
}

/**
 * @brief Calculates the second-order gradient (curvature) of the polynomial function at a specific point.
 *
 * @param regression_state Pointer to the RegressionState structure containing the coefficients.
 * @param x The point at which to calculate the second-order gradient.
 * @return The second-order gradient at the given point x.
 */
double calculate_second_order_gradient(const RegressionState *regression_state, double x)
{
    DEBUG_PRINT_3("Calculating second-order gradient at x=%.6f\n", x);

    double second_derivative = 0.0;
    double x_power = 1.0;
    int degree = regression_state->polynomial_degree;

    for (int i = 2; i <= degree; ++i)
    {
        x_power *= x;
        second_derivative += i * (i - 1) * regression_state->coefficients[i] * x_power / (x * x);
        DEBUG_PRINT_3("Term %d: i=%d, coefficient=%.6f, x_power=%.6f, second_derivative=%.6f\n", i, i, regression_state->coefficients[i], x_power, second_derivative);
    }

    DEBUG_PRINT_3("Second-order gradient calculated: %.6f\n", second_derivative);
    return second_derivative;
}

/**
 * @brief Tracks the values (phase angles) added to the RLS array and calculates gradients after each addition.
 *
 * This generalized function calculates gradients (first-order, second-order, etc.) for a given dataset
 * using Recursive Least Squares (RLS) polynomial regression. It accepts a function pointer to determine
 * the type of gradient to calculate, making it flexible for various gradient computations.
 *
 * @param measurements Array of data points, each of type MqsRawDataPoint_t. We use the phaseAngle from each element.
 * @param length The number of points to add starting from the given start index.
 * @param start_index The index in the measurements array from which to begin adding values.
 * @param degree The degree of the polynomial to use for regression (e.g., 2 for quadratic, 3 for cubic).
 * @param calculate_gradient Function pointer to the gradient calculation function (e.g., first-order, second-order).
 * @param result Pointer to store the gradient calculation results.
 */
void trackGradients(
    const MqsRawDataPoint_t *measurements,
    uint16_t length,
    uint16_t start_index,
    uint8_t degree,
    double (*calculate_gradient)(const RegressionState *, double),
    GradientCalculationResult *result)
{
    DEBUG_PRINT_3("Entering trackGradients with startIndex=%u, length=%u, degree=%u\n",
                  start_index, length, degree);

    // Initialize the result structure
    result->size = 0;
    DEBUG_PRINT_2("Initialized GradientCalculationResult: size=0\n");

    // Initialize the regression state
    RegressionState regression_state;
    initialize_regression_state(&regression_state, degree, RLS_WINDOW);
    DEBUG_PRINT_3("RegressionState initialized with degree=%u\n", degree);

    // Loop through the measurements and calculate gradients
    for (uint16_t i = 0; i < length; ++i)
    {
        uint16_t current_index = start_index + i;

        // Extract the phase angle from the MqsRawDataPoint_t struct
        double measurement = (double)measurements[current_index].phaseAngle;

        // Add the phase angle to the regression
        // (Depending on your existing API for 'add_data_point_to_regression',
        //  you may need to adapt it to accept a single double instead of an array.)
        add_data_point_to_regression(&regression_state, measurements, current_index);
        DEBUG_PRINT_1("Added measurement %.6f to RegressionState (index=%u)\n",
                      measurement, current_index);

        // Get the x_value corresponding to the current data point
        double x_value = (double)(regression_state.total_data_points_added - 1);
        DEBUG_PRINT_2("Current x_value: %.2f (total_data_points_added=%u)\n",
                      x_value, regression_state.total_data_points_added);

        // Check if there are enough points to perform regression and calculate gradient
        if (regression_state.current_num_points >= (regression_state.polynomial_degree + 1))
        {
            // Calculate the gradient using the provided function pointer
            double gradient = calculate_gradient(&regression_state, x_value);
            DEBUG_PRINT_1("Calculated gradient: %.6f at x=%.2f\n", gradient, x_value);

            // Store the gradient in the result struct if there's space
            if (result->size < RLS_WINDOW)
            {
                result->gradients[result->size++] = gradient;
                DEBUG_PRINT_2("Stored gradient in result->gradients[%u]\n", result->size - 1);
            }
            else
            {
                // Handle overflow: log a warning and stop collecting gradients
                DEBUG_PRINT_1("Gradient array overflow at index %u. Maximum window size reached.\n", result->size);
                break;
            }
        }
        else
        {
            // Not enough points to calculate gradient
            DEBUG_PRINT_1("Insufficient points to calculate gradient (current_num_points=%u)\n",
                          regression_state.current_num_points);
        }
    }

    // Post-processing after collecting gradients
    if (result->size > 0)
    {
        DEBUG_PRINT_2("Collected %u gradients\n", result->size);
    }
    else
    {
        DEBUG_PRINT_2("No gradients were collected\n");
    }

    // Note: Median and MAD calculations are handled in the calling function (e.g., identifyTrends)
    DEBUG_PRINT_3("Exiting trackGradients\n");
}

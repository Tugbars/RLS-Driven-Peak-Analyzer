#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "windowed_running_gradient.h"
#include "running_cubic_gradient.h"
#include "running_peak_analysis.h"
#include "running_quadratic_gradient.h"

#include "rls_analysis_parameters.h"

int main() {
	const double data[] = { 11.26, 11.13, 11.276, 11.136, 11.194, 11.22, 11.18, 11.371, 11.269, 11.389, 11.377, 11.448, 11.248, 11.569, 11.442, 11.332, 11.525, 11.532, 11.427, 11.48, 11.431, 11.435, 11.515, 11.561, 11.489, 11.575, 11.693, 11.484, 11.619, 11.572, 11.708, 11.676, 11.528, 11.828, 11.598, 11.697, 11.823, 11.839, 11.889, 11.839, 11.927, 11.99, 11.826, 11.964, 11.813, 12.089, 12.076, 11.986, 12.192, 12.118, 11.997, 12.144, 12.311, 12.315, 12.293, 12.441, 12.35, 12.354, 12.399, 12.419, 12.455, 12.453, 12.566, 12.538, 12.725, 12.588, 12.89, 12.914, 12.913, 12.94, 13.023, 12.842, 13.141, 13.083, 13.321, 13.352, 13.372, 13.488, 13.555, 13.674, 13.66, 13.651, 13.997, 13.868, 13.969, 14.105, 14.191, 14.531, 14.278, 14.432, 14.523, 14.744, 14.441, 14.802, 15.066, 14.889, 15.088, 15.493, 15.494, 15.774, 16.081, 15.847, 15.77, 16.463, 16.188, 16.436, 16.388, 16.708, 16.842, 16.95, 17.413, 17.781, 17.929, 18.59, 18.478, 18.806, 18.738, 19.276, 19.744, 19.83, 20.157, 20.391, 20.367, 20.866, 21.19, 21.535, 21.994, 22.578, 22.869, 22.826, 23.654, 23.704, 24.724, 24.866, 25.609, 25.97, 26.53, 26.882, 28.051, 28.477, 29.039, 29.041, 30.509, 30.461, 31.413, 31.605, 32.46, 34.372, 34.084, 35.04, 35.014, 35.321, 36.64, 36.967, 37.323, 38.137, 39.721, 39.566, 40.457, 41.6, 42.555, 42.642, 43.222, 43.819, 43.611, 44.176, 44.248, 44.378, 44.385, 44.624, 44.413, 44.228, 43.884, 43.298, 43.342, 42.24, 41.29, 41.042, 40.263, 39.176, 38.162, 38.101, 36.628, 37.511, 35.735, 34.998, 33.44, 32.395, 31.241, 30.892, 29.863, 29.021, 28.65, 27.799, 26.718, 27.519, 25.865, 25.085, 23.494, 23.679, 22.736, 22.747, 22.706, 21.48, 21.577, 21.855, 20.304, 20.415, 19.763, 19.755, 19.377, 18.424, 18.938, 18.784, 18.147, 18.072, 17.786, 17.284, 17.527, 17.015, 16.833, 16.676, 16.693, 16.288, 16.071, 15.955, 15.876, 15.482, 15.065, 15.203, 15.227, 14.815, 15.08, 14.961, 14.7, 14.835, 14.526, 14.142, 14.214, 14.261, 14.06, 14.094, 13.824, 13.916, 13.689, 13.736, 13.817, 13.489, 13.69, 13.645, 13.332, 13.344, 13.16, 13.13, 13.111, 12.94, 12.997, 12.829, 12.906, 12.569, 12.824, 12.533, 12.456, 12.47, 12.578, 12.363, 12.328, 12.399, 12.285, 12.3, 12.28, 12.232, 12.366, 12.303, 12.174, 12.035, 11.957, 12.123, 11.931, 12.031, 11.943, 12.024, 11.989, 11.93, 11.702, 11.943, 11.827, 11.818, 11.877, 11.696, 11.784, 11.726, 11.617, 11.542, 11.503, 11.58, 11.639, 11.688, 11.514, 11.541, 11.267, 11.388, 11.45, 11.537, 11.489, 11.32, 11.366, 11.376, 11.269, 11.24, 11.417, 11.314, 11.296, 11.306, 11.231, 11.381, 11.173, 11.321, 11.21, 11.185, 11.298, 11.121, 11.287, 11.227, 11.112, 11.199, 11.208, 11.224, 11.21, 11.168, 11.168, 11.266, 11.075, 11.15, 10.992, 11.005, 11.081, 10.916, 10.984, 11.074, 10.954, 11.052, 11.105, 10.999, 10.953, 11.02, 10.945, 11.056, 11.065, 10.913, 10.958, 11.022, 11.038, 10.969, 10.887, 10.904, 10.936, 10.899, 10.988, 10.752 };
	size_t size = sizeof(data) / sizeof(data[0]);
	/*
	 for (int i = 0; i < size; i++) {
	 printf("%d. %f ", i, data[i]);
	 if ((i + 1) % 10 == 0) {
	     printf("\n");
	 }
	}
	*/
	
	init_cubic_rls_analysis_parameters(
        7.0,   // significance_thresh: Threshold for determining significant cubic trends.
               // Trends with a sum of gradients above this threshold are considered significant.
        5,     // duration_thresh: Minimum number of consecutive points required for a trend to be considered significant.
               // Ensures that only sustained trends are analyzed.
        2,     // min_trend_count: Minimum number of consistent trends required for cubic analysis.
               // Helps filter out noise and minor fluctuations.
        2,     // max_third_order_trend_decrease_count: Maximum number of consecutive decreases allowed in a cubic trend.
               // Allows for minor fluctuations without discarding a potentially valid trend.
        2      // max_third_order_trend_increase_count: Maximum number of consecutive increases allowed in a cubic trend.
               // Helps to filter out noise when identifying a cubic trend.
    );
    
    // Initialize the quadratic RLS analysis parameters with the necessary thresholds
    init_quadratic_rls_analysis_parameters(
        0.3,   // min_gradient_sum: Minimum cumulative sum of second-order gradients required for a trend to be considered valid.
               // Helps in identifying significant concave or convex regions.
        2,     // max_decrease_count: Maximum number of consecutive negative gradients allowed when tracking an increasing trend.
               // Allows for minor fluctuations without discarding a potentially valid trend.
        2      // max_increase_count: Maximum number of consecutive positive gradients allowed when tracking a decreasing trend.
               // Helps to filter out noise when identifying a decreasing trend.
    );
    
    // Initialize the on-peak analysis parameters for average gradient thresholds and consistent trend count
    init_on_peak_analysis_parameters(
        0.45,  // min_avg_increase: Minimum average increase required to consider a significant increasing trend.
               // Ensures that only trends with substantial average gradient are flagged as significant.
        -0.45, // min_avg_decrease: Minimum average decrease required to consider a significant decreasing trend.
               // Ensures that only trends with substantial average gradient are flagged as significant.
        5      // min_consistent_trend_count: Minimum number of consecutive data points required for a trend to be considered consistent.
               // Helps ensure that only sustained trends are analyzed.
    );

	segment_trend_and_concavity_analysis(data, size, 140, 0.5);

	return 0;
}

//Segment analysis çok iyi değil. üzerinde durmam lazım. 
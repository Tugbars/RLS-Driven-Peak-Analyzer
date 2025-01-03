/**
 * @file sliding_window_analysis_overview.h
 * @brief High-level overview and rationale behind the sliding window analysis for resonant peak detection.
 *
 * @details
 * This documentation provides an outline of the **adaptive sliding window approach** used to find
 * and precisely center on a resonant peak in a frequency-domain dataset. Instead of performing
 * one massive sweep over the entire frequency range, the algorithm focuses on a smaller analysis
 * window, which it slides or expands incrementally based on measured gradients (from RLS polynomial
 * regression). By doing so, the program adaptively homes in on the likely peak region and verifies
 * if the peak is genuinely centered. 
 *
 * The key elements of this strategy include:
 * - **Recursive Least Squares (RLS) polynomial regression** to smooth out noise and obtain robust
 *   gradients (first- and second-order).
 * - **Gradient-based decisions** about whether the window indicates an upward slope, a downward slope,
 *   or a peak region.
 * - A **state machine** that orchestrates window movement, peak centering, and final validation.
 *
 * Below is a comprehensive summary of the approach:
 *
 * ### 1. Overview
 * - **Goal**: Find and precisely center on a resonant peak in a frequency-domain dataset without doing
 *   a single massive sweep.
 * - **Approach**:
 *   1. **Collect a small window of data** (phase angles at specific frequencies).
 *   2. **Use RLS** to smooth/fit a polynomial to these points, computing gradients (slopes, curvature).
 *   3. **Analyze** whether the window is rising, falling, or contains a peak.
 *   4. **Shift the window** left/right (or expand it) based on gradient analysis.
 *   5. **Verify** if the peak is properly centered (strong rise on the left, strong fall on the right).
 *   6. **Repeat** until a confirmed, centered peak is found (or until boundary/time limits are reached).
 *
 * ### 2. Why a Sliding Window?
 * - **Data Can Be Noisy**:
 *   - Large sweeps can be cumbersome. Focusing on a smaller slice ("window") at a time allows for more
 *     efficient polynomial regression. It also keeps the analysis "local," which is useful for honing in
 *     on the actual peak.
 * - **Adaptive Sweeping**:
 *   - By not scanning blindly over every frequency, the algorithm can quickly home in on where the peak
 *     likely is, potentially saving time and handling unknown damping factors or broad/narrow peaks more gracefully.
 *
 * ### 3. RLS Polynomial Regression
 * - **Incremental Updates**:
 *   - New data points are added (old points removed) without recalculating the entire polynomial fit from
 *     scratch, thanks to Recursive Least Squares.
 * - **Rolling Polynomial Fit**:
 *   - Maintains a "rolling" polynomial fit of the current window. Once fit, it calculates:
 *     - **First-order gradient** (slope) to see if we are going up or down.
 *     - **Second-order gradient** (curvature) to detect the shape of the potential peak.
 *
 * ### 4. Determining Peak vs. No Peak
 * - **Gradient Threshold**:
 *   - The algorithm checks for gradients above a specific threshold to determine slope significance, 
 *     accounting for variations in material damping.
 * - **Centering the Peak**:
 *   - If a peak is detected but not centered, the window is shifted based on gradient trends. 
 *   - A peak is considered centered if consistent negative gradients lead to a value of approximately -1.0.
 *   - The centering process continues until the peak is fully captured and satisfies these conditions.
 * - **Peak Verification**:
 *   - Once the peak is suspected, the algorithm verifies it by checking the **second-order gradients**:
 *     - The left side of the window should show strongly positive second-order gradients.
 *     - The right side of the window should show strongly negative second-order gradients.
 * - **Handling Truncation**:
 *   - If the analysis suggests part of the peak is outside the current window, the window is expanded or shifted
 *     to fully encompass the peak.
 *
 * ### 5. Error Control
 * To ensure the algorithm doesn’t get stuck in oscillations:
 * - If the sliding window moves left, then right, and left again, or moves right, then left, and right again, 
 *   the analysis is **cancelled**, as this indicates the algorithm is stuck in an ambiguous state.
 * - Once an error condition is detected, the process resets or terminates based on program rules.
 *
 * ### 6. Buffer Management
 * - The impedance analyzer chip outputs values that are stored in a **local buffer**.
 * - This buffer begins filling from the **middle** of the buffer array and dynamically expands:
 *   - To the left if the sliding window moves toward lower frequencies.
 *   - To the right if the sliding window moves toward higher frequencies.
 * - This design ensures efficient memory usage while maintaining flexibility to adapt to analysis requirements.
 *

 * ### 7 DirectionDecisions How the Program Decides to Move Left or Right
 * The function `determineMoveDirection()` decides the sliding window's movement based on gradient trends:
 * 1. **Trend Counts**: Analyzes the longest consistent increasing and decreasing trends to assess data distribution.
 * 2. **Threshold Comparison**: Moves if trend counts exceed `TREND_THRESHOLD`, signaling significant slopes.
 * 3. **Global Gradients Check**: Evaluates maximum and minimum gradients to detect if values fall within an "undecided" range.
 * 4. **Final Decision**:
 *    - Strong increase only: Move right.
 *    - Strong decrease only: Move left.
 *    - Both strong: Compare gradient sums to determine the dominant direction.
 *    - Weak or balanced: Remain `UNDECIDED`.
 * 
 * ### 7. State Machine Logic
 * The entire adaptive process is implemented as a **state machine**:
 * - **WAITING**: Idle until a new sweep starts.
 * - **INITIAL_ANALYSIS**: Load the initial window, perform a quick RLS fit, and analyze initial gradients.
 * - **SEGMENT_ANALYSIS / UPDATE_BUFFER_DIRECTION**: Determine the slope direction and shift the window if necessary.
 * - **UNDECIDED_TREND_CASE**: Handle ambiguous data patterns by making small adjustments and re-checking.
 * - **PEAK_CENTERING**: Attempt to center the peak if it’s partially visible within the window.
 * - **PEAK_FINDING_ANALYSIS**: Verify if the peak satisfies all gradient-based conditions.
 * - **PEAK_TRUNCATION_HANDLING**: Expand or shift the window if the peak is partially outside its bounds.
 * - **EXPAND_ANALYSIS_WINDOW**: Increase the number of data points in the window if more resolution is needed.
 * - **WAITING (again)**: Return to idle once the analysis is complete or if error limits are reached.
 *
 * @note
 * This adaptive process ensures robust, noise-resistant detection of resonant peaks using gradient-based 
 * decisions and dynamic window adjustments, all orchestrated through state machine transitions.
 *
 * ### 9. In Short
 * - **Small Windows, Smoothed**: RLS polynomial fits ensure stable gradient estimates within small data slices.
 * - **Adaptive**: The algorithm dynamically shifts or expands the window to find the peak efficiently.
 * - **Verification**: Additional checks (second-order gradients, consistent slopes) confirm peak validity and centering.
 * - **Error Handling**: Mechanisms detect and reset the process if the algorithm becomes stuck.
 * - **State Machine Control**: A clear series of states governs the analysis flow, ensuring systematic decisions.
 *
 * @author
 * Tugbars Heptaskin
 * @date
 * 2024-12-30
 *
 * @version
 * 1.0
 */


#include "sliding_window_analysis.h"
#include "ics.h"
//#include <hal/ics/ics.h>
#include <stdio.h>


#define FORGETTING_FACTOR_ADJUSTMENT 0.2

#define CENTERING_RATIO 2

#define PEAK_VERIFICATION_COUNT 5

#define MAX_SAVGOL_WINDOW_SIZE 25

/******************************************************************************/
/* Global Variables */
/******************************************************************************/

/**
 * @brief Union representing the status of the sliding window analysis process.
 *
 * This union contains flags that indicate the current status of the analysis,
 * such as whether a peak is found, if the analysis is undecided, or if the peak is centered.
 */
typedef union {
	uint8_t value; /**< The combined value of all status flags. */
	struct {
		uint8_t isPeakFound : 1;            /**< Flag indicating if a peak is found. */
		uint8_t isUndecided : 1;            /**< Flag indicating if the analysis is undecided. */
		uint8_t isCentered : 1;             /**< Flag indicating if the peak is centered. */
		uint8_t isSweepRequested : 1;       /**< Flag indicating if a sweep is requested. */
		uint8_t isSweepDone : 1;            /**< Flag indicating if the sweep is done. */
		uint8_t isNotCentered : 1;          /**< Flag indicating if the peak is still not centered. */
		uint8_t isVerificationFailed : 1;   /**< Flag indicating if peak verification failed. */
		uint8_t isBoundaryError : 1;        /**< Flag indicating if a boundary error occurred. */
	};
} SwpStatus_t;

/** @brief Context for sliding window analysis. */
SlidingWindowAnalysisContext ctx; // Sliding window analysis context

/** @brief Current status of the sliding window analysis. */
SwpStatus_t currentStatus = { .value = 0x0 };  // Initialize all flags to 0

/** @brief Current state of the sliding window analysis state machine. */
SwpState_t currentState = SWP_WAITING;  // Initialize to SWP_WAITING state

// Add a global flag to track boundary error
bool boundaryErrorOccurred = false;

/*******************************************************************************
 * Undecided deadlock error control
 ******************************************************************************/

/**
 * @brief Tracks the number of consecutive times the state machine enters the undecided case.
 *un
 * This variable is incremented each time the state machine enters the `SWP_UNDECIDED_TREND_CASE` state.
 * It limits the number of times the system can enter the undecided state within a session. If the counter
 * exceeds a predefined limit (e.g., 7), the system will transition to the `SWP_WAITING` state to prevent
 * the state machine from getting stuck in an undecided loop. This helps ensure that the analysis progresses
 * or halts gracefully if a decision cannot be reached.
 */
static uint8_t undecided_case_counter = 0;

/**
 * @brief Flag indicating that an error occurred due to exceeding the undecided case limit.
 *
 * This flag is set to `true` when the `undecided_case_counter` exceeds the maximum allowed value
 * (e.g., 7 undecided states). When this flag is set, the system will transition to the `SWP_WAITING` state,
 * and a warning will be displayed in the callback function executed by the state machine. The flag
 * will be reset when the `SWP_WAITING` state is entered, indicating that the system is ready for
 * a new session or sweep.
 */
bool undecided_error_flag = false;

/** @brief Maximum number of times the state machine can enter the undecided case within a session.
 *
 * This variable defines the limit on how many times the state machine can enter the `SWP_UNDECIDED_TREND_CASE`
 * state before it triggers an error. If the `undecided_case_counter` exceeds this value, the system will
 * transition to the `SWP_WAITING` state and raise an error flag, indicating that a resolution could not be reached.
 */
static const uint8_t MAX_UNDECIDED_CASE_ATTEMPTS = 7;

/*******************************************************************************
 * Peak centering deadlock error control
 ******************************************************************************/

/**
 * @brief Defines the maximum number of peak centering attempts.
 *
 * This constant sets a limit on the number of times the algorithm will attempt to center the peak.
 * After reaching this limit, the system will stop further attempts to prevent an infinite loop in cases
 * where the peak cannot be accurately centered. This ensures that the process either succeeds within
 * a reasonable number of attempts or gracefully exits if centering is not possible.
 */
#define MAX_CENTERING_ATTEMPTS 4

/**
 * @brief Tracks the number of attempts made to center the peak.
 *
 * This variable is used to limit the number of centering attempts in the peak centering state.
 * Each time the peak centering process fails to accurately center the peak, this counter is
 * incremented. Once it reaches a predefined maximum number of attempts (MAX_CENTERING_ATTEMPTS),
 * the system will exit the centering state and transition to a waiting or error state to prevent
 * an infinite loop in cases where the peak cannot be centered properly.
 */
static uint8_t centering_attempts = 0;

/**
 * @brief Adjusts the forgetting factor used in peak centering analysis.
 *
 * The forgetting factor controls how much weight is given to more recent data points during the
 * gradient analysis for peak centering. This variable is incremented by 0.1 after each centering
 * attempt to give more weight to the newer data points, allowing the algorithm to make progressively
 * larger adjustments in subsequent attempts. This helps improve the accuracy of peak centering
 * over multiple attempts.
 */
// static double centering_forgetting_factor = 0.7;

/*******************************************************************************
 * Repetitive shift deadlock error control
 ******************************************************************************/

/**
* @brief Size of the shift tracker array used to detect alternating left and right shifts.
*
* This constant defines the size of the `shift_tracker` array, which stores the last two
* shift directions (LEFT, RIGHT) made during the sliding window analysis. The purpose is to
* detect whether the analysis is getting stuck in an alternating pattern of shifts.
*/
#define SHIFT_TRACKER_SIZE 2

/**
 * @brief Array to track the last two shift directions.
 *
 * This array stores the most recent shift directions (LEFT or RIGHT) that occurred during the
 * sliding window analysis. If the shift directions alternate between LEFT and RIGHT over the
 * last two shifts, the system concludes that the analysis is stuck and raises an error.
 */
static PeakPosition shift_tracker[SHIFT_TRACKER_SIZE] = {UNDECIDED, UNDECIDED};

/**
 * @brief Index for tracking the current position in the `shift_tracker` array.
 *
 * This variable tracks the current index in the `shift_tracker` array, which stores the last two
 * shift directions. The array functions as a circular buffer, where the index wraps around after
 * reaching the size of the array (`SHIFT_TRACKER_SIZE`).
 */
static uint8_t shift_tracker_index = 0;

/**
 * @brief Updates the shift tracker with the latest direction.
 *
 * This function updates the `shift_tracker` array with the most recent shift direction. The array
 * functions as a circular buffer, meaning that once it reaches its size limit, it starts overwriting
 * the oldest direction. The function also checks whether the last two shift directions alternated
 * between LEFT and RIGHT. If this condition is met, an error is raised, and the analysis transitions
 * to the `SWP_WAITING` state.
 *
 * @param new_direction The new shift direction (LEFT, RIGHT, or UNDECIDED).
 */
static void update_shift_tracker(PeakPosition new_direction) {
	// Update the shift tracker with the latest direction
	shift_tracker[shift_tracker_index] = new_direction;
	shift_tracker_index = (shift_tracker_index + 1) % SHIFT_TRACKER_SIZE;

	// Check if the analysis is stuck (LEFT b RIGHT b LEFT or RIGHT b LEFT b RIGHT)
	if ((shift_tracker[0] == LEFT_SIDE && shift_tracker[1] == RIGHT_SIDE) ||
	        (shift_tracker[0] == RIGHT_SIDE && shift_tracker[1] == LEFT_SIDE)) {
		printf("Error: Sliding window analysis is stuck alternating between left and right shifts.\n");
		boundaryErrorOccurred = true;  // Raise an error flag
		currentStatus.isSweepDone = 1; // End the analysis and go to SWP_WAITING
	}
}

/******************************************************************************/
/* Function Prototypes (Internal Functions) */
/******************************************************************************/

/* OnEntry and OnExit function declarations for each state */
static void OnEntryInitialAnalysis(void);
static void OnEntrySegmentAnalysis(void);
static void OnEntryUpdateBufferDirection(void);
static void OnEntryUndecidedTrendCase(void);
static void OnEntryPeakFindingAnalysis(void);
static void OnEntryWaiting(void);
static void OnEntryPeakTruncationHandling(void);
static void OnEntryPeakCentering(void);
static void OnEntryExpandAnalysisWindow(void);

static void OnExitInitialAnalysis(void);
static void OnExitSegmentAnalysis(void);
static void OnExitUpdateBufferDirection(void);
static void OnExitUndecidedTrendCase(void);
static void OnExitPeakFindingAnalysis(void);
static void OnExitWaiting(void);
static void OnExitPeakCentering(void);
static void OnExitPeakTruncationHandling(void);
static void OnExitExpandAnalysisWindow(void);

// Forward declaration of startAdaptiveSweep function
static void startAdaptiveSweep(void);
static void adaptiveSweepSampleCb(double real, double imaginary, bool isSweepComplete); //static void adaptiveSweepSampleCb(int16_t real, int16_t imaginary, bool isSweepComplete);
void SwpProcessStateChange(void);
static void initBufferManager(MqsRawDataPoint_t* dataBuffer, int start_index);
static SwpState_t NextState(SwpState_t state);

/******************************************************************************/
/* State Function Definitions */
/******************************************************************************/

/**
 * @brief Structure containing function pointers and completion status for state functions.
 */
typedef struct {
	void (*onEntry)(void);
	void (*onExit)(void);
	bool isComplete;
} StateFuncs_t;

/**
 * @brief Array of state functions and their completion status.
 *
 * This array maps each state to its corresponding onEntry and onExit functions, along with a flag
 * indicating whether the state's processing is complete.
 */
static StateFuncs_t STATE_FUNCS[SWP_STATE_LAST] = {
	{OnEntryInitialAnalysis, OnExitInitialAnalysis, false},         // SWP_INITIAL_ANALYSIS
	{OnEntrySegmentAnalysis, OnExitSegmentAnalysis, false},         // SWP_SEGMENT_ANALYSIS
	{OnEntryUpdateBufferDirection, OnExitUpdateBufferDirection, false}, // SWP_UPDATE_BUFFER_DIRECTION
	{OnEntryUndecidedTrendCase, OnExitUndecidedTrendCase, false},   // SWP_UNDECIDED_TREND_CASE
	{OnEntryPeakCentering, OnExitPeakCentering, false},             // SWP_PEAK_CENTERING
	{OnEntryPeakFindingAnalysis, OnExitPeakFindingAnalysis, false}, // SWP_PEAK_FINDING_ANALYSIS
	{OnEntryPeakTruncationHandling, OnExitPeakTruncationHandling, false}, // SWP_PEAK_TRUNCATION_HANDLING
	{OnEntryExpandAnalysisWindow, OnExitExpandAnalysisWindow, false}, // SWP_EXPAND_ANALYSIS_WINDOW
	{OnEntryWaiting, OnExitWaiting, false}                          // SWP_WAITING
};
/******************************************************************************/
/* Sliding Window Analysis Functions */
/******************************************************************************/

/**
 * @brief This function concludes the sweep process by determining whether to continue
 * with adaptive sweep or proceed with the next state change. If buffer updates are needed,
 * it calls `startAdaptiveSweep`. Otherwise, it proceeds to call `SwpProcessStateChange`.
 */
static void concludeSweepState(void) {
	if (buffer_update_info.needs_update) {
		startAdaptiveSweep();  // Take note: This is the adaptive sweep function that will be called
	} else {

		SwpProcessStateChange();
	}
}

/**
 * @brief This function is the callback for handling sweep data collection during an adaptive sweep process.
 *
 * It simulates reading data from the AD5933 chip and processes the real and imaginary values
 * by passing them to `AdptSweepAddDataPoint`.
 *
 * @param real The real part of the complex data (simulated as a double).
 * @param imaginary The imaginary part of the complex data (simulated as a double).
 * @param isSweepComplete Flag indicating whether the sweep process is complete.
 */
static void adaptiveSweepSampleCb(double real, double imaginary, bool isSweepComplete) {
	// Simulate data collection by calling `AdptSweepAddDataPoint` with the new real and imaginary values.
	AdptSweepAddDataPoint(real, imaginary);  // Pass the real and imaginary values directly to the function

	// Reset the update flag after processing
	buffer_update_info.needs_update = false;

	// After completing data collection, finalize the sweep process
	if (isSweepComplete) {
		currentStatus.isSweepDone = true;
		// Trigger the state change after the sweep is done
		SwpProcessStateChange();
	}
}


/**
 * @brief This function initiates an adaptive sweep process, activating the AD5933 chip and collecting data.
 * After collecting data, it calls `concludeSweepState` to either continue or finalize the sweep.
 */
static void startAdaptiveSweep(void) {
	//set it to 0 before starting a sweep. the sweep count is used in mapping the collected data points with the data buffer.
	currentRawSweep->count = 0;

	HalIcsSetStartFreq(buffer_update_info.phase_index_start);
	HalIcsSetFreqInc(currentRawSweep->setup->frequencyIncrement);
	HalIcsSetNoOfSamples(buffer_update_info.move_amount);
	HalIcsStartFreqSweep(adaptiveSweepSampleCb);
}


/**
 * @brief Initializes the buffer manager for sliding window analysis.
 *
 * This function sets up the buffer manager with the appropriate parameters, including the buffer size,
 * window size, starting index, and frequency parameters.
 */
static void initBufferManager(MqsRawDataPoint_t* dataBuffer, int start_index) {
	// Initialize the buffer manager with appropriate parameters
	// Arguments:
	// buffer          -> The buffer that will store the phase and impedance values (array of MqsRawDataPoint_t)
	// BUFFER_SIZE     -> Total size of the buffer (maximum number of data points the buffer can hold)
	// 30              -> Sliding window size (number of data points to analyze in each window)
	// 0               -> Starting index in the phaseAngles array (this could be any index where you want to start the analysis)
	// 11300.0         -> Starting frequency for the frequency sweep (in Hz, e.g., starting at 11300 Hz)
	// 1.0             -> Frequency increment per step (in Hz, e.g., increment by 1 Hz for each data point)
	init_buffer_manager(dataBuffer, MQS_SWEEP_MAX_NUMBER_OF_SAMPLES, RLS_WINDOW, start_index, 11300.0, 1.5);

	// Debugging: Print initialization state
	//printf("[DEBUG] Buffer Manager initialized:\n");
	//printf("Buffer size: %d, Window size: %d, Start index: %d\n", BUFFER_SIZE, WINDOW_SIZE, start_index);
}


/**
 * @brief Starts the sliding window analysis process.
 *
 * This function initializes the context for the sliding window analysis, sets up the buffer manager,
 * and starts the state machine for processing.
 *
 * @param phaseAngles Pointer to the array of phase angles.
 * @param phase_angle_size Size of the phase angle array.
 * @param callback Function pointer to the callback function to be executed after analysis.
 */
void startSlidingWindowAnalysis(MesSweep_t *sweep, int start_index, Callback_t callback) { //ILK MES START OLARAK GOREV ALABILIR.
	ctx.callback = callback;
	ctx.isTruncatedLeft = false;
	ctx.isTruncatedRight = false;
	ctx.peakIndex = 0;

	// Initialize the buffer manager
	initBufferManager(sweep->data, start_index);

	// Enable sweep request
	currentStatus.isSweepRequested = 1;  // Start sweep request

	// Set the initial state to waiting
	currentState = SWP_WAITING;

	// Set the initial state and start the state machine
	SwpProcessStateChange();  // Start the state machine
}


/******************************************************************************/
/* State Machine Entry Functions */
/******************************************************************************/

/**
 * @brief Entry function for the SWP_INITIAL_ANALYSIS state.
 *
 * This function loads the initial buffer with phase angle data and sets the sweep request flag.
 */
static void OnEntryInitialAnalysis(void) {
	// Call load_initial_buffer to set up buffer update
	load_initial_buffer();

	// Set the sweep request flag and mark this state as complete
	currentStatus.isSweepRequested = true;
	STATE_FUNCS[SWP_INITIAL_ANALYSIS].isComplete = true;

	// Transition to conclude the sweep state, which will handle the async sweep
	concludeSweepState();
}


/**
 * @brief Entry function for the SWP_SEGMENT_ANALYSIS state.
 *
 * This function performs segment analysis on the current window of data to determine the direction of movement.
 * If the analysis is stuck alternating between left and right shifts, it raises an error and transitions to
 * the `SWP_WAITING` state to prevent infinite looping.
 */
static void OnEntrySegmentAnalysis(void) {

	uint8_t degree = 3;             // Degree for RLS regression
	GradientOrder gradientOrder = GRADIENT_ORDER_FIRST;  // Change as needed

	//printf("=== Performing Peak Analysis ===\n");
	PeakAnalysisResult peakResult = performPeakAnalysis(
	                                    buffer_manager.buffer,
	                                    RLS_WINDOW, //not necessary
	                                    buffer_manager.current_buffer_index,
	                                    buffer_manager.window_size,
	                                    degree,
	                                    gradientOrder
	                                );

	// isSweepDone can be activated if the peak is verified in performPeakAnalysis. For debugging purposes I do not see the point.
	// the state machine is however built around it. 

	// peakResult.isCenteringNeeded, OnPeak utilize edilmiyor sorun o.
	if(peakResult.isSignificantPeak)
	{
		currentStatus.isPeakFound = 1;
	}

	// Update the context with the direction result
	ctx.direction = peakResult.moveDirection;

	// Update the shift tracker with the new direction
	update_shift_tracker(ctx.direction);

	// Check if we are on the peak
	if (ctx.direction == ON_THE_PEAK) { 
	    
	} else if (ctx.direction == UNDECIDED || ctx.direction == NEGATIVE_UNDECIDED) {
		currentStatus.isUndecided = 1;  // Raise undecided flag for both cases
	}

	// Mark this state as complete and proceed to the next state
	STATE_FUNCS[SWP_SEGMENT_ANALYSIS].isComplete = true;
	SwpProcessStateChange();
}


/**
 * @brief Entry function for the SWP_UPDATE_BUFFER_DIRECTION state.
 *
 * This function updates the buffer based on the determined direction from the segment analysis.
 */
static void OnEntryUpdateBufferDirection(void) {

	// Set up buffer update info
	update_buffer_for_direction(ctx.direction, buffer_manager.window_size / 2);

	STATE_FUNCS[SWP_UPDATE_BUFFER_DIRECTION].isComplete = true;

	// After setting up, perform the buffer update if needed
	concludeSweepState();
}


/**
 * @brief Entry function for the SWP_UNDECIDED_TREND_CASE state.
 *
 * This function handles the case when the segment analysis result is undecided by moving the window forward.
 * If it's a negative undecided case, it moves the window backward instead.
 */
static void OnEntryUndecidedTrendCase(void) {

	if (undecided_case_counter >= MAX_UNDECIDED_CASE_ATTEMPTS) {
		// Set the error flag and transition to SWP_WAITING
		printf("Undecided case limit exceeded. Setting error flag and returning to SWP_WAITING.\n");
		undecided_error_flag = true;
		currentState = SWP_WAITING;
		SwpProcessStateChange();
		return;
	}

	if (ctx.direction == NEGATIVE_UNDECIDED) {
		handle_negative_undecided_case();  // Move the window backward
	} else {
		handle_undecided_case();  // Move the window forward
	}

	STATE_FUNCS[SWP_UNDECIDED_TREND_CASE].isComplete = true;

	// After setting up, perform the buffer update if needed
	concludeSweepState();
}


/**
 * @brief Resets state completion flags and checks if the maximum centering attempts have been exceeded.
 *
 * This function resets the necessary flags for the peak centering state to ensure that
 * the state is not prematurely marked as complete. It also checks if the number of
 * centering attempts has reached the predefined maximum (`MAX_CENTERING_ATTEMPTS`). If so,
 * the sweep is marked as complete, and the state machine transitions to the next state.
 *
 * ### Steps:
 * 1. Resets the `isComplete` flag for the current state to indicate the processing has started.
 * 2. Resets the `isNotCentered` flag for tracking the centering status.
 * 3. Checks if the number of centering attempts has exceeded `MAX_CENTERING_ATTEMPTS`.
 *    - If the maximum is reached, it marks the sweep as done and moves the state machine
 *      to the next phase.
 *
 * @note
 * - If the peak cannot be centered within the maximum number of attempts, the state transitions
 *   to SWP_WAITING and marks the sweep as done to prevent further iterations.
 *
 * @see SwpProcessStateChange
 */
static void reset_centering_flags(void) {
	STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = false;
	currentStatus.isNotCentered = 0;

	if (centering_attempts >= MAX_CENTERING_ATTEMPTS) {
		printf("Max centering attempts reached. Transitioning to WAITING state.\n");
		currentStatus.isSweepDone = 1;
		STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
		SwpProcessStateChange();
	}
}


/**
 * @brief Converts a MoveDirection enum to the corresponding integer argument
 *        for move_window_and_update_if_needed.
 *
 * @param moveDirection The MoveDirection enum value to convert.
 * @return int The corresponding integer argument for the move_window_and_update_if_needed function.
 */
static int convertMoveDirectionToArgument(MoveDirection moveDirection) {
	switch (moveDirection) {
	case MOVE_LEFT:
		return LEFT_SIDE;  // Assuming LEFT_SIDE is defined as the integer expected
	case MOVE_RIGHT:
		return RIGHT_SIDE; // Assuming RIGHT_SIDE is defined as the integer expected
	default:
		return -1;  // Invalid direction; handle this case appropriately
	}
}

/**
 * @brief Entry function for the SWP_PEAK_CENTERING state, orchestrating the peak centering process.
 *
 * This function serves as the main entry point for peak centering within the SWP_PEAK_CENTERING state.
 * It coordinates the workflow by calling various helper functions that manage different aspects of
 * the centering process, including resetting state flags, performing gradient calculations, analyzing
 * trends, and adjusting the data window.
 *
 * ### Function Workflow:
 * 1. **Reset Flags**: Calls `reset_centering_flags` to ensure that the state is properly initialized and checks
 *    whether the maximum number of centering attempts has been exceeded.
 * 2. **Perform Peak Centering Logic**: Computes the total second-order gradient to analyze the data curvature.
 * 3. **Check Centering Status**: Verifies if the peak is centered by comparing the gradient sum to a threshold.
 *    - If centered, transitions to the next state using `SwpProcessStateChange`.
 * 4. **Analyze Trends**: Uses quadratic regression to track increasing and decreasing trends.
 * 5. **Handle Invalid Trend Data**: If trend data is invalid, terminates the centering process early.
 * 6. **Adjust Window Position**: If the peak is not centered, adjusts the data window based on the trend analysis.
 * 7. **Finalize Attempt**: Finalizes the centering attempt by marking the state as complete and performing any remaining tasks.
 *
 * @note
 * - The function orchestrates the flow by invoking other helper functions that perform specific tasks such as gradient calculation, trend analysis, and window adjustment.
 * - The actual logic for handling each step (e.g., computing gradients, adjusting the window) is handled by the external functions.
 *
 * @see reset_centering_flags
 * @see perform_peak_centering_logic
 * @see check_centering_status
 * @see track_gradient_trends_with_quadratic_regression
 * @see handle_invalid_trend_data
 * @see adjust_window_position
 * @see finalize_centering_attempt
 */
static void OnEntryPeakCentering(void) {
	// Reset flags and check if centering attempts exceeded
	reset_centering_flags();
	if (STATE_FUNCS[SWP_PEAK_CENTERING].isComplete) return; //where doe sthis come from?


	uint8_t degree = 3;             // Degree for RLS regression
	printf("\n=== Centering ===\n");
	// Call identifyAndCalculateCentering to get centering parameters
	PeakCenteringResult centeringResult = identifyAndCalculateCentering(
	        buffer_manager.buffer,
	        buffer_manager.current_buffer_index,
	        buffer_manager.buffer_size,
	        degree
	                                      );
	// Handle invalid trend data
	//if (handle_invalid_trend_data(&gradient_trends)) return;

	// Convert moveDirection to the appropriate argument for move_window_and_update_if_needed
	int directionArg = convertMoveDirectionToArgument(centeringResult.moveDirection);

	if(!centeringResult.isCentered) {
		move_window_and_update_if_needed(directionArg,  centeringResult.moveAmount);
		//case SWP_PEAK_FINDING_ANALYSIS: will determine if this is really centered or not, after this function.
	}

	currentStatus.isCentered = 1;
	STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
	concludeSweepState();  // Finalize buffer updates if needed
}

/**
 * @brief Handles and logs warnings based on peak truncation flags.
 *
 * This function evaluates the truncation flags in the given `QuadraticPeakAnalysisResult` structure
 * and logs warnings about left or right truncation. Additionally, it updates a provided context
 * structure (e.g., `ctx`) with truncation information if necessary.
 *
 * @param result Pointer to the `QuadraticPeakAnalysisResult` containing verification results.
 * @param ctx    Pointer to a user-defined context structure to store truncation state (optional).
 */
void handle_truncation_warnings(const QuadraticPeakAnalysisResult *result, SlidingWindowAnalysisContext *ctx) {
	if (!result) {
		printf("Invalid result pointer provided to handle_truncation_warnings.\n");
		return;
	}

	// Handle left truncation warning
	if (result->is_truncated_left) {
		printf("Warning: Peak verification truncated on the left side.\n");
		if (ctx) {
			ctx->isTruncatedLeft = true;  // Update context if provided
		}
	}

	// Handle right truncation warning
	if (result->is_truncated_right) {
		printf("Warning: Peak verification truncated on the right side.\n");
		if (ctx) {
			ctx->isTruncatedRight = true;  // Update context if provided
		}
	}

	// If no truncations occurred
	if (!result->is_truncated_left && !result->is_truncated_right) {
		printf("No truncation detected in peak verification.\n");
	}
}


/**
 * @brief Processes the result of peak verification and determines if the peak is valid.
 *
 * This function evaluates the result of a peak verification process. It resets relevant counters,
 * updates the status flags, and marks the sweep as complete if the peak is valid. If the peak
 * verification fails, the function flags that the peak is not centered and prompts a return to
 * the peak-centering process.
 *
 * ### Key Operations:
 * - Resets centering attempts upon successful peak verification.
 * - Marks the sweep as complete if the peak is valid.
 * - Flags the peak as not centered if verification fails.
 * - Logs detailed messages for success or failure scenarios.
 *
 * @param result Pointer to the `QuadraticPeakAnalysisResult` containing verification results.
 * @param centering_attempts Pointer to the centering attempts counter to reset if successful.
 * @return `true` if the peak was successfully verified, `false` otherwise.
 * @see QuadraticPeakAnalysisResult
 */
/**
 * @brief Processes the result of peak verification and determines if the peak is valid.
 *
 * This function evaluates the result of a peak verification process. It resets relevant counters,
 * updates the status flags, and marks the sweep as complete if the peak is valid. If the peak
 * verification fails, the function flags that the peak is not centered and prompts a return to
 * the peak-centering process.
 *
 * ### Key Operations:
 * - Resets centering attempts upon successful peak verification.
 * - Marks the sweep as complete if the peak is valid.
 * - Flags the peak as not centered if verification fails.
 * - Logs detailed messages for success or failure scenarios.
 *
 * @param result Pointer to the `QuadraticPeakAnalysisResult` containing verification results.
 * @return `true` if the peak was successfully verified, `false` otherwise.
 * @see QuadraticPeakAnalysisResult
 */
bool process_peak_verification(
    const QuadraticPeakAnalysisResult *result
) {
	if (!result) {
		printf("Error: Null pointer passed to process_peak_verification.\n");
		return false;
	}

	// If peak verification succeeded
	if (result->peak_found) {
		printf("[Success] Peak verification successful. Peak is centered at index %u.\n", result->peak_index);

		// Reset counters and mark sweep as done
		centering_attempts = 0;                  // Reset the centering attempts counter
		currentStatus.isSweepDone = 1;           // Mark the sweep as done
		currentStatus.isNotCentered = 0;         // Clear the not-centered flag
		currentStatus.isVerificationFailed = 0;   // Clear the verification failed flag

		return true;
	} else {
		// If peak verification failed
		printf("[Failure] Peak verification failed. Returning to peak centering.\n");

		// Update status flags
		currentStatus.isNotCentered = 1;         // Set flag to indicate the peak is not centered
		currentStatus.isVerificationFailed = 1;  // Set flag to indicate verification failure

		return false;
	}
}

/**
 * @brief Executes the peak verification process using quadratic regression.
 *
 * This function verifies the peak by performing a quadratic regression analysis over the
 * current buffer data. It checks for the presence of a valid peak and handles any truncation
 * scenarios. If the peak is valid, it prints the analysis interval and proceeds with the
 * state transition. If not, it flags the peak as not centered.
 *
 * ### Key Operations:
 * - Calls `find_and_verify_quadratic_peak` to analyze the peak.
 * - Prints truncation warnings if the analysis is truncated.
 * - Processes the peak verification result.
 *
 * @see find_and_verify_quadratic_peak
 * @see print_truncation_warnings
 * @see process_peak_verification
 */
static uint16_t perform_peak_verification(void) {
	// Perform peak verification through quadratic regression
	int degree = 3;

	QuadraticPeakAnalysisResult verificationResult = verifyPeakValidity(
	            buffer_manager.buffer,
	            buffer_manager.buffer_size,
	            buffer_manager.current_buffer_index,
	            RLS_WINDOW,  //another problem.
	            degree
	        );

	//printf("peak index of the peak: %u\n", verificationResult.peak_index);

	//Print truncation warnings, if any
	handle_truncation_warnings(&verificationResult, &ctx); // handle_truncation_warnings

	// Process peak verification result, if successful or not
	if (process_peak_verification(&verificationResult)) {  //process_peak_verification
        ctx.peakIndex = verificationResult.peak_index;
		return verificationResult.peak_index; // Peak is centered, no further action needed
	}
}


/**
 * @brief Checks if the peak is centered and performs verification if centered.
 *
 * This function coordinates the entire process of verifying if the peak is centered.
 * It first calculates the total second-order gradient sum to check the peak's centration.
 * If the peak is centered, the function proceeds to perform the peak verification process.
 *
 * ### Key Operations:
 * - Computes the total sum of second-order gradients to assess centration.
 * - Verifies the peak using quadratic regression if the peak is centered.
 *
 * @see compute_total_second_order_gradient
 * @see perform_peak_verification
 * @see verify_peak_centering
 */
static void check_peak_centering_and_verify(void) { //trackGradients would solve it.
	// Compute total second-order gradient sum to verify peak

	// Initialize result structure
	GradientCalculationResult result;

	// Define minimum ratio (e.g., 0.3 for 30%)
	double min_ratio = 0.3;

	bool isCentered = peakAnalysis_checkPeakCentering(  // daha dikkat etmek gerekiyor buna.
	                      buffer_manager.buffer,
	                      buffer_manager.buffer_size,
	                      buffer_manager.current_buffer_index,
	                      min_ratio,
	                      &result
	                  );

	if (isCentered) {
		printf("Valid negative trend with adequate peak capture detected.\n");
		printf("Start Index: %u, End Index: %u\n", result.startIndex, result.endIndex);
		currentStatus.isCentered = 1;
		perform_peak_verification(); // should return the index.
		
	} else {
		printf("Negative trend not found or peak not adequately captured.\n");
		currentStatus.isCentered = 0;
	}
}


/**
 * @brief Entry function for the SWP_PEAK_FINDING_ANALYSIS state.
 *
 * This function performs several key tasks to verify the centering and validity of the detected peak:
 *
 * ### Step 1: Re-check the Centration of the Peak
 * - The function first computes the total sum of the second-order gradients across the current data window
 *   using a quadratic regression analysis. This sum reflects the curvature of the data and helps determine
 *   whether the peak is centered in the window.
 * - If the sum is outside the predefined threshold range (`centered_gradient_sum`), it indicates that the peak is not yet
 *   centered, and the function returns the state machine to the peak-centering state.
 *
 * ### Step 2: Perform Peak Verification
 * - If the sum of second-order gradients is within the acceptable threshold, the function proceeds to verify the peak.
 * - The verification process, done via `find_and_verify_quadratic_peak`, checks if the data exhibits a valid peak shape:
 *   - Increasing trends on the left side of the peak.
 *   - Decreasing trends on the right side of the peak.
 * - It also accounts for truncation scenarios. If the verification process encounters the boundaries of the window before
 *   completing the necessary trend checks, truncation flags (`is_truncated_left` and `is_truncated_right`) are set.
 *
 * ### Step 3: Handle Truncation Scenarios
 * - If the peak verification is truncated on either the left or right side, this information is logged.
 * - The function prints specific messages indicating whether truncation occurred on the left, right, or both sides
 *   of the data window.
 *
 * ### Step 4: Final Peak Verification and State Transition
 * - If the peak is successfully verified based on the trend analysis, the function prints the sliding window analysis
 *   interval and proceeds to mark the sweep as complete (`isSweepDone` flag).
 * - If the peak verification fails (e.g., the data does not exhibit the required trend patterns), the function returns
 *   to the peak-centering state by setting the `isNotCentered` flag.
 *
 * ### Summary of Key Functionality:
 * - Re-checks peak centration based on the sum of second-order gradients.
 * - Verifies the peak using quadratic regression and trend analysis.
 * - Handles boundary truncations during verification.
 * - Directs the state machine to either the peak-centering or peak-verification complete state.
 *
 * @see compute_total_second_order_gradient
 * @see find_and_verify_quadratic_peak
 * @see print_analysis_interval
 * @see SwpProcessStateChange
 */
static void OnEntryPeakFindingAnalysis(void) {
	check_peak_centering_and_verify();

	// Mark the state as complete and transition to the next state
	STATE_FUNCS[SWP_PEAK_FINDING_ANALYSIS].isComplete = true;
	SwpProcessStateChange();  // Process the next state
}


/**
 * @brief Entry function for the SWP_PEAK_TRUNCATION_HANDLING state.
 *
 * This function handles the scenario where peak verification during the peak finding analysis
 * is incomplete due to data truncation on the left or right side of the data window.
 * It adjusts the data window by moving it in the appropriate direction and attempts to
 * re-verify the peak after the adjustment.
 *
 * ### Function Workflow:
 * 1. **Initialization**:
 *    - Initializes a flag `peak_verified_after_truncation` to track if the peak is verified
 *      after handling truncation.
 * 2. **Left Side Truncation Handling**:
 *    - Checks if the peak verification was truncated on the left side (`ctx.isTruncatedLeft`).
 *    - If so, it moves the data window to the left by a predefined amount (e.g., 5 positions)
 *      using `move_window_and_update_if_needed`.BI
 *    - Re-verifies the peak at the current buffer index by calling `verify_peak_at_index`.
 *    - Resets the `isTruncatedLeft` flag in the context.
 * 3. **Right Side Truncation Handling**:
 *    - Checks if the peak verification was truncated on the right side (`ctx.isTruncatedRight`).
 *    - If so, it moves the data window to the right by a predefined amount (e.g., 5 positions)
 *      using `move_window_and_update_if_needed`.
 *    - Re-verifies the peak at the current buffer index by calling `verify_peak_at_index`.
 *    - Resets the `isTruncatedRight` flag in the context.
 * 4. **Peak Verification Post-Truncation**:
 *    - If the peak is successfully verified after adjustments, it marks the sweep as done
 *      by setting `currentStatus.isSweepDone`.
 *    - If the peak verification fails after adjustments, it sets `currentStatus.isVerificationFailed`
 *      to indicate that the peak verification failed, prompting a return to peak finding analysis.
 * 5. **State Transition**:
 *    - Marks the current state as complete (`isComplete = true`).
 *    - Calls `SwpProcessStateChange()` to transition to the next state based on the updated status flags.
 *
 * ### Related Functions:
 * - **`find_and_verify_quadratic_peak`**:
 *   - Detects and verifies peaks using second-order gradients from quadratic regression.
 *   - Returns a `QuadraticPeakAnalysisResult` containing peak verification status and truncation flags.
 * - **`verify_quadratic_peak`**:
 *   - Verifies a detected peak by checking for consistent increasing trends on the left side and
 *     decreasing trends on the right side of the peak.
 *   - Considers truncation scenarios if the data window does not fully encompass the required trends.
 *
 * ### Notes:
 * - The function relies on the context (`ctx`) to check for truncation flags (`isTruncatedLeft`, `isTruncatedRight`)
 *   and to reset them after handling.
 * - The function modifies `currentStatus` to communicate the result of the truncation handling to the state machine.
 * - The data window adjustments are crucial for ensuring that the peak is fully captured within the analysis window,
 *   which is essential for accurate peak verification.
 */
static void OnEntryPeakTruncationHandling(void) {
	// Since this state may involve multiple transitions based on truncation handling,
	// we keep isComplete = false initially
	STATE_FUNCS[SWP_PEAK_TRUNCATION_HANDLING].isComplete = false;

	if (ctx.isTruncatedLeft) {
		printf("Handling truncation on the left side.\n");
		move_window_and_update_if_needed(LEFT_SIDE, PEAK_VERIFICATION_COUNT);

		// Perform buffer update if needed
		concludeSweepState();
	}

	if (ctx.isTruncatedRight) {
		printf("Handling truncation on the right side.\n");
		move_window_and_update_if_needed(RIGHT_SIDE, PEAK_VERIFICATION_COUNT);

		// Perform buffer update if needed
		concludeSweepState();
	}

	// Don't mark the state complete here since we will evaluate the outcome in OnExit
}

/**
 * @brief Entry function for the SWP_EXPAND_ANALYSIS_WINDOW state.
 *
 * This state ensures enough data points are available for the Savitzky-Golay filter.
 * It adjusts the buffer window by expanding it left or right as needed.
 */
static void OnEntryExpandAnalysisWindow(void) {
	//printf("[SWP_EXPAND_ANALYSIS_WINDOW] Expanding the analysis window for Savitzky-Golay filter.\n");

	// Check proximity to boundaries using the new function
	BoundaryProximityResult proximityResult = checkBoundaryProximity(ctx.peakIndex, MAX_SAVGOL_WINDOW_SIZE / 2);

	// Handle proximity to boundaries
	if (proximityResult.isNearBoundary) {
		printf("[SWP_EXPAND_ANALYSIS_WINDOW] Boundary proximity detected! Direction: %s, Move Amount: %d\n",
		       proximityResult.direction == LEFT_SIDE ? "LEFT" : "RIGHT",
		       proximityResult.moveAmount);

		// Expand the window based on the detected proximity
		if (proximityResult.direction == LEFT_SIDE) {
			move_window_and_update_if_needed(LEFT_SIDE, proximityResult.moveAmount);
		} else if (proximityResult.direction == RIGHT_SIDE) {
			move_window_and_update_if_needed(RIGHT_SIDE, proximityResult.moveAmount);
		}
	} else {
		//printf("[SWP_EXPAND_ANALYSIS_WINDOW] No boundary proximity detected. No window expansion needed.\n");
	}
	
	print_analysis_interval();

	// Mark the state as complete
	STATE_FUNCS[SWP_EXPAND_ANALYSIS_WINDOW].isComplete = true;

	// Transition to SWP_WAITING after expanding the window
	concludeSweepState();
}

/**
 * @brief Entry function for the SWP_WAITING state.
 *
 * This function waits for a sweep request and executes the callback if provided.
 */
static void OnEntryWaiting(void) {

	if (currentStatus.isSweepRequested) {
		STATE_FUNCS[SWP_WAITING].isComplete = true;
	}
}

/******************************************************************************/
/* Exit Functions */
/******************************************************************************/
static void OnExitInitialAnalysis(void) {
	//printf("Exiting SWP_INITIAL_ANALYSIS state.\n");
}

static void OnExitSegmentAnalysis(void) {
	//printf("Exiting SWP_SEGMENT_ANALYSIS state.\n");
}

static void OnExitUpdateBufferDirection(void) {
	//printf("Exiting SWP_UPDATE_BUFFER_DIRECTION state.\n");
}

static void OnExitUndecidedTrendCase(void) {
	//printf("Exiting SWP_UNDECIDED_TREND_CASE state.\n");
}

static void OnExitPeakFindingAnalysis(void) {
	//printf("Exiting SWP_PEAK_FINDING_ANALYSIS state.\n");
}

static void OnExitWaiting(void) {
	//printf("Exiting SWP_WAITING state.\n");
}

static void OnExitPeakCentering(void) {
	//printf("Exiting SWP_PEAK_CENTERING state.\n");
}

/**
 * @brief Exit function for the SWP_PEAK_TRUNCATION_HANDLING state.
 *
 * This function handles the peak verification after truncation and buffer update are done.
 */
static void OnExitPeakTruncationHandling(void) {
	// Verify if the peak was successfully centered after truncation

	// Now call verifyPeakValidity with the new startIndex
	QuadraticPeakAnalysisResult verificationResult = verifyPeakValidity(
	            buffer_manager.buffer,
	            buffer_manager.buffer_size,
	            buffer_manager.current_buffer_index,
	            RLS_WINDOW,  //another problem.
	            3 //degree
	        );


	if (verificationResult.peak_found && !verificationResult.is_truncated_left && !verificationResult.is_truncated_right) {
		printf("Peak verification successful after truncation handling, peak is centered.\n");
	

		ctx.peakIndex = verificationResult.peak_index;
		currentStatus.isSweepDone = 1;  // Mark the sweep as done
		// Now the state is ready to transition, mark it complete

	} else {
		// If peak is not verified, return to peak finding analysis
		printf("Peak verification failed after truncation handling, returning to peak finding analysis.\n");
		currentStatus.isVerificationFailed = 1;  // Set flag to indicate verification failed
	}

	STATE_FUNCS[SWP_PEAK_TRUNCATION_HANDLING].isComplete = true;
}

/**
 * @brief Exit function for the SWP_EXPAND_ANALYSIS_WINDOW state.
 */
static void OnExitExpandAnalysisWindow(void) {
	//printf("[SWP_EXPAND_ANALYSIS_WINDOW] Exiting state.\n");
}

/******************************************************************************/
/* State Machine Process and Transitions */
/******************************************************************************/

/**
 * @brief Processes state changes in the sliding window analysis state machine.
 *
 * This function manages the state transitions in the sliding window analysis state machine.
 * It ensures that the appropriate OnExit and OnEntry functions are called during state transitions.
 * The state machine continues processing until no further state changes occur.
 *
 * ### Structure and Purpose:
 * - The function uses a loop to repeatedly check if a state transition is needed.
 * - Each state is represented by a set of functions: `OnEntry`, `OnExit`, and a `isComplete` flag that indicates
 *   whether the state has finished its processing.
 * - State transitions are triggered based on the `NextState` function, which determines the appropriate state to move to
 *   after completing the current one.
 *
 * ### Key Concepts:
 * - **State Completion (`isComplete`)**: Each state has an `isComplete` flag to indicate whether the state's processing
 *   is done. The state machine will only transition to the next state if this flag is set. If the state is not marked
 *   as complete, the next state is often determined by the logic inside the `OnExit()` function of the current state,
 *   allowing flexibility for more complex state transitions that are dependent on the exit logic of a state.
 *
 * - **OnExit/OnEntry Functions**:
 *   - When a state is marked as incomplete (`!isComplete`), the `OnExit()` function (if defined) is called to handle
 *     cleanup, and potentially dictate what the next state will be.
 *   - Once a state transition is made, the `OnEntry()` function for the new state is invoked to initialize the new state's processing.
 *
 * - **State Transition (`NextState`)**: The `NextState()` function determines the next state to transition into
 *   based on the current state and flags. This function, in combination with the `OnExit()` function, ensures the
 *   correct state is selected, especially for states that require further processing or specific conditions before moving forward.
 *
 * ### Intention:
 * - The function aims to provide flexibility in state transitions by allowing states to either fully complete and move on
 *   to the next state, or re-enter themselves if they require further processing, such as in asynchronous operations.
 * - The `OnExit()` logic allows states to determine the appropriate next state during exit, particularly useful in scenarios
 *   where a state's exit depends on external factors or complex logic.
 * - The loop ensures the state machine processes transitions in sequence without breaking the flow, making it adaptable for real-time or multi-step operations.
 *
 * ### Working Mechanism:
 * 1. **State Completion Check**:
 *    - If the current state is incomplete (`!STATE_FUNCS[currentState].isComplete`), the function calls the `OnExit()` function (if defined)
 *      to handle cleanup and determine any next steps.
 *    - The `NextState()` function is called to determine the next state, and `stateChanged` tracks if the state has changed.
 *    - This allows a state to finalize any pending tasks in its `OnExit()` before transitioning to the next one.
 *
 * 2. **State Transition**:
 *    - After calling `OnExit()`, the `isComplete` flag for the new state is reset, and the `OnEntry()` function for the new state is invoked.
 *    - The `currentStatus.value` is reset to ensure no status flags carry over to the next state.
 *
 * 3. **Complete State Handling**:
 *    - If a state is marked as complete (`STATE_FUNCS[currentState].isComplete`), the function checks for the next state and processes it.
 *    - The loop continues processing state transitions until no more state changes occur.
 *
 * 4. **Callback Handling**:
 *    - If the state machine reaches the `SWP_WAITING` state and a callback is set (`ctx.callback`), the callback is executed and cleared.
 *
 * ### Example Use Case:
 * The function is suitable for scenarios where:
 * - Some states need multiple passes to complete their processing (e.g., waiting for an asynchronous event).
 * - Transitions between states depend on conditions evaluated during the `OnExit()` function.
 * - A callback function is executed after the state machine reaches a specific terminal state (e.g., `SWP_WAITING`).
 *
 * ### Flexibility:
 * The function handles both synchronous and asynchronous transitions. States that require multiple iterations
 * to complete can use the `OnExit()` logic to determine the next state dynamically. Meanwhile, states that do not
 * require re-entry can transition immediately.
 *
 * @note The state machine continues processing until no further state changes occur, ensuring all transitions are properly handled.
 *
 */
void SwpProcessStateChange(void) {
	bool stateChanged;

	// The state machine should process until no further state changes occur
	do {
		// If the current state is incomplete, handle its exit
		if (!STATE_FUNCS[currentState].isComplete) {
			// Ensure this debug message is only printed during truncation handling or incomplete states

			// Execute OnExit if it's not complete
			if (STATE_FUNCS[currentState].onExit != NULL) {
				STATE_FUNCS[currentState].onExit();
			}

			// Determine the next state after OnExit
			SwpState_t nextState = NextState(currentState);
			stateChanged = (nextState != currentState);

			// Reset the currentStatus before transitioning
			currentStatus.value = 0x0;

			if (stateChanged) {

				// Ensure the old state's completion flag is set to false
				STATE_FUNCS[currentState].isComplete = false;
				// Move to the new state
				currentState = nextState;

				STATE_FUNCS[currentState].isComplete = false;

				// Call the new state's OnEntry function
				if (STATE_FUNCS[currentState].onEntry != NULL) {
					STATE_FUNCS[currentState].onEntry();
				}
			}
		}

		// Process the state if it is marked complete
		if (STATE_FUNCS[currentState].isComplete) {
			//printf("[DEBUG] State %d is complete. Proceeding to next state...\n", currentState);

			// Determine the next state
			SwpState_t nextState = NextState(currentState);
			stateChanged = (nextState != currentState);

			// Clear currentStatus flags
			currentStatus.value = 0x0;

			if (stateChanged) {
				// Transition to the new state
				currentState = nextState;

				// Ensure the new state's completion flag is set to false
				STATE_FUNCS[currentState].isComplete = false;

				// Call the OnEntry function for the new state
				if (STATE_FUNCS[currentState].onEntry != NULL) {
					STATE_FUNCS[currentState].onEntry();
				}
			}
		}

	} while (stateChanged && STATE_FUNCS[currentState].isComplete);

	// Handle the callback execution when we reach the waiting state
	if (currentState == SWP_WAITING && ctx.callback) {
		printf("[DEBUG] Executing callback as we have reached SWP_WAITING.\n");
		ctx.callback();  // Execute the callback function
		ctx.callback = NULL;  // Clear the callback to avoid re-execution
	}
}

/******************************************************************************/
/* Helper Functions for State Transitioning */
/******************************************************************************/

/**
 * @brief Determines the next state in the state machine based on the current state and status flags.
 *
 * This function is used to manage state transitions in the sliding window analysis process. It checks for various flags
 * such as boundary errors, sweep completion, and centering conditions to determine the appropriate next state.
 *
 * @param state The current state.
 * @return SwpState_t The next state to transition to.
 */
static SwpState_t NextState(SwpState_t state) {
	switch (state) {
	case SWP_INITIAL_ANALYSIS:
		if (currentStatus.isPeakFound) {
			return SWP_PEAK_CENTERING;
		}
		if (currentStatus.isUndecided) {
			return SWP_UNDECIDED_TREND_CASE;
		}
		return SWP_SEGMENT_ANALYSIS;

	case SWP_SEGMENT_ANALYSIS:
		if (currentStatus.isUndecided) {
			return SWP_UNDECIDED_TREND_CASE;
		}
		if (currentStatus.isPeakFound) {
			return SWP_PEAK_CENTERING;
		}
		if (currentStatus.isSweepDone) {
			return SWP_WAITING;  // Transition to WAITING state if sweep is done
		}
		return SWP_UPDATE_BUFFER_DIRECTION;

	case SWP_UPDATE_BUFFER_DIRECTION:
		if (currentStatus.isBoundaryError) {
			return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
		}
		return SWP_SEGMENT_ANALYSIS;

	case SWP_PEAK_CENTERING:
		if (currentStatus.isBoundaryError) {
			return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
		}
		if (currentStatus.isCentered) {
			return SWP_PEAK_FINDING_ANALYSIS;
		}
		if (currentStatus.isSweepDone) {
			return SWP_EXPAND_ANALYSIS_WINDOW;
		}
		return SWP_PEAK_CENTERING;

	case SWP_UNDECIDED_TREND_CASE:
		// Check if undecided case counter has reached the limit
		if (undecided_case_counter >= MAX_UNDECIDED_CASE_ATTEMPTS) {
			return SWP_WAITING;  // Transition to SWP_WAITING if limit exceeded
		}
		if (currentStatus.isBoundaryError) {
			return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
		}
		return SWP_SEGMENT_ANALYSIS;

	case SWP_PEAK_FINDING_ANALYSIS:
		if (currentStatus.isBoundaryError) {
			return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
		}
		if (ctx.isTruncatedLeft || ctx.isTruncatedRight) {
			return SWP_PEAK_TRUNCATION_HANDLING;
		}
		if (!currentStatus.isCentered) {
			return SWP_PEAK_CENTERING;
		}
		if (currentStatus.isSweepDone) {
			return SWP_EXPAND_ANALYSIS_WINDOW;
		}
		return SWP_PEAK_FINDING_ANALYSIS;

	case SWP_PEAK_TRUNCATION_HANDLING:
		if (currentStatus.isBoundaryError) {
			return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
		}
		if (currentStatus.isSweepDone) {
			return SWP_EXPAND_ANALYSIS_WINDOW; // new
		}
		if (currentStatus.isVerificationFailed) {
			return SWP_PEAK_FINDING_ANALYSIS;
		}
		return SWP_PEAK_FINDING_ANALYSIS;
		
	case SWP_EXPAND_ANALYSIS_WINDOW:
            if (STATE_FUNCS[SWP_EXPAND_ANALYSIS_WINDOW].isComplete) {
                return SWP_WAITING;  // Transition to SWP_WAITING when this state is complete
            }
            return SWP_EXPAND_ANALYSIS_WINDOW;

	case SWP_WAITING:
		if (currentStatus.isSweepRequested) {
			return SWP_INITIAL_ANALYSIS;
		}
		return SWP_WAITING;

	default:
		return SWP_WAITING;
	}
}

/**
 * @brief Sets the boundary error flag in the currentStatus.
 *
 * This function allows external modules (like buffer_manager) to set the boundary error flag in the sliding window
 * analysis state machine, which will trigger a state transition to SWP_WAITING.
 *
 * @param flag The value to set for the boundary error flag (1 for error, 0 to clear).
 */
// Function to set the boundary error flag
void set_boundary_error_flag(uint8_t flag) {
	currentStatus.isBoundaryError = flag;
	if (flag) {
		boundaryErrorOccurred = true; // Set the persistent boundary error flag
	}
}


/**
 * @brief Computes and returns the analysis interval information for Savitzky-Golay filtering.
 *
 * This function combines the adjusted buffer index and the absolute distance
 * between the analysis start and end indices into a struct.
 *
 * @return SavgolWindowInterval A struct containing the adjusted buffer index and
 *         the analysis interval distance. If indices are invalid, returns -1 for both fields.
 */
SavgolWindowInterval get_savgol_window_interval(void) {
    SavgolWindowInterval interval;

    // Calculate the adjusted buffer index
    if (analysis_start_index == -1 || analysis_end_index == -1) {
        interval.adjustedBufferIndex = -1;
        interval.analysisIntervalDistance = -1;
    } else {
        int relative_position = analysis_start_index - buffer_manager.current_phase_index;
        interval.adjustedBufferIndex = (buffer_manager.current_buffer_index + relative_position - buffer_shift_offset + buffer_manager.buffer_size) % buffer_manager.buffer_size;

        // Calculate the absolute distance between analysis indices
        interval.analysisIntervalDistance = abs(analysis_end_index - analysis_start_index);
    }

    return interval;
}


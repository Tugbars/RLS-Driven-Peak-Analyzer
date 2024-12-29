# RLS based Adaptive Sweep for a handheld Impedance Analyzer

## **High-Level Overview**

This documentation provides an outline of the **adaptive sliding window approach** used to find and precisely center on a resonant peak in a frequency-domain dataset. Instead of performing one massive sweep over the entire frequency range, the algorithm focuses on a smaller analysis window, which it slides or expands incrementally based on measured gradients (from RLS polynomial regression). By doing so, the program adaptively homes in on the likely peak region and verifies if the peak is genuinely centered.

### **Key Elements**
- **Recursive Least Squares (RLS) polynomial regression** to smooth out noise and obtain robust gradients (first- and second-order).
- **Gradient-based decisions** about whether the window indicates an upward slope, a downward slope, or a peak region.
- A **state machine** that orchestrates window movement, peak centering, and final validation.

---

## **1. Overview**
- **Goal**: Find and precisely center on a resonant peak in a frequency-domain dataset without doing a single massive sweep.
- **Approach**:
  1. **Collect a small window of data** (phase angles at specific frequencies).
  2. **Use RLS** to smooth/fit a polynomial to these points, computing gradients (slopes, curvature).
  3. **Analyze** whether the window is rising, falling, or contains a peak.
  4. **Shift the window** left/right (or expand it) based on gradient analysis.
  5. **Verify** if the peak is properly centered (strong rise on the left, strong fall on the right).
  6. **Repeat** until a confirmed, centered peak is found (or until boundary/time limits are reached).

---

## **2. Why a Sliding Window?**
### **Data Can Be Noisy**
- Large sweeps can be cumbersome. Focusing on a smaller slice ("window") at a time allows for more efficient polynomial regression. It also keeps the analysis "local," which is useful for honing in on the actual peak.

### **Adaptive Sweeping**
- By not scanning blindly over every frequency, the algorithm can quickly home in on where the peak likely is, potentially saving time and handling unknown damping factors or broad/narrow peaks more gracefully.

---

## **3. RLS Polynomial Regression**
- **Incremental Updates**:
  - New data points are added (old points removed) without recalculating the entire polynomial fit from scratch, thanks to Recursive Least Squares.
- **Rolling Polynomial Fit**:
  - Maintains a "rolling" polynomial fit of the current window. Once fit, it calculates:
    - **First-order gradient** (slope) to see if we are going up or down.
    - **Second-order gradient** (curvature) to detect the shape of the potential peak.

---

## **4. Determining Peak vs. No Peak**
- **Gradient Threshold**:
  - The algorithm checks for gradients above a specific threshold to determine slope significance, accounting for variations in material damping.
- **Centering the Peak**:
  - If a peak is detected but not centered, the window is shifted based on gradient trends.
  - A peak is considered centered if consistent negative gradients lead to a value of approximately `-1.0`.
- **Peak Verification**:
  - Verifies the peak by checking the **second-order gradients**:
    - Left side: strongly positive second-order gradients.
    - Right side: strongly negative second-order gradients.
- **Handling Truncation**:
  - If part of the peak is outside the window, the window is expanded or shifted to fully encompass it.

---

## **5. Error Control**
To ensure the algorithm doesn’t get stuck:
- If the sliding window moves left, then right, and left again, or moves right, then left, and right again, the analysis is **cancelled**.
- Once an error condition is detected, the process resets or terminates based on program rules.

---

## **6. Buffer Management**
- The impedance analyzer chip outputs values that are stored in a **local buffer**.
- This buffer begins filling from the **middle** of the buffer array and dynamically expands:
  - To the left if the sliding window moves toward lower frequencies.
  - To the right if the sliding window moves toward higher frequencies.
- This design ensures efficient memory usage while maintaining flexibility.

---

## **7. How the Program Decides to Move Left or Right**
The function `determineMoveDirection()` decides the sliding window's movement based on gradient trends:

1. **Trend Counts**:
   - Analyzes the longest consistent increasing and decreasing trends to assess data distribution.
2. **Threshold Comparison**:
   - Moves if trend counts exceed `TREND_THRESHOLD`, signaling significant slopes.
3. **Global Gradients Check**:
   - Evaluates maximum and minimum gradients to detect if values fall within an "undecided" range.
4. **Final Decision**:
   - **Strong increase only**: Move right.
   - **Strong decrease only**: Move left.
   - **Both strong**: Compare gradient sums to determine the dominant direction.
   - **Weak or balanced**: Remain `UNDECIDED`.

---

## **8. State Machine Logic**
The entire adaptive process is implemented as a **state machine**:
- **WAITING**: Idle until a new sweep starts.
- **INITIAL_ANALYSIS**: Load the initial window, perform a quick RLS fit, and analyze initial gradients.
- **SEGMENT_ANALYSIS / UPDATE_BUFFER_DIRECTION**: Determine the slope direction and shift the window if necessary.
- **UNDECIDED_TREND_CASE**: Handle ambiguous data patterns by making small adjustments and re-checking.
- **PEAK_CENTERING**: Attempt to center the peak if it’s partially visible within the window.
- **PEAK_FINDING_ANALYSIS**: Verify if the peak satisfies all gradient-based conditions.
- **PEAK_TRUNCATION_HANDLING**: Expand or shift the window if the peak is partially outside its bounds.
- **EXPAND_ANALYSIS_WINDOW**: Increase the number of data points in the window if more resolution is needed.
- **WAITING (again)**: Return to idle once the analysis is complete or if error limits are reached.
![ActivityUML](https://github.com/user-attachments/assets/e128460d-4f67-44c8-ac2f-869e4a8a98f8)

---

## **9. Summary**
- **Small Windows, Smoothed**: RLS polynomial fits ensure stable gradient estimates within small data slices.
- **Adaptive**: The algorithm dynamically shifts or expands the window to find the peak efficiently.
- **Verification**: Additional checks (second-order gradients, consistent slopes) confirm peak validity and centering.
- **Error Handling**: Mechanisms detect and reset the process if the algorithm becomes stuck.
- **State Machine Control**: A clear series of states governs the analysis flow, ensuring systematic decisions.

---

### Author
**Tugbars Heptaskin**  
Date: **2024-12-30**  
Version: **1.0**

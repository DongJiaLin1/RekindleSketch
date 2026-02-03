#ifndef PARA_H
#define PARA_H

#include <string>

/**
 * Configuration parameters for RekindleSketch
 * Adjust these values to configure memory usage and algorithm parameters
 */

// ===== Memory Configuration =====
const size_t TOTAL_MEMORY_BYTES = 200 * 1024;  // Total memory budget (200KB)
const int DEFAULT_N = 550;                      // Number of cells per bucket
const size_t ESTIMATED_CELL_SIZE = 16;          // Estimated memory per cell (bytes)

// Calculate number of buckets based on memory constraints
const int CALCULATED_M = static_cast<int>(
    TOTAL_MEMORY_BYTES / (ESTIMATED_CELL_SIZE * DEFAULT_N)
);

const int DEFAULT_M = (CALCULATED_M > 0) ? CALCULATED_M : 1;  // Number of buckets

// ===== Memory Score Parameters =====
const double DEFAULT_ALPHA = 0.015;  // Exponential decay rate
const double DEFAULT_BETA = 1.5;     // Reward function parameter
const double DEFAULT_Q = 1.0;        // Base reward value
const int DEFAULT_Z = 1;             // Replacement probability parameter

// ===== Sliding Window Parameters =====
const int R = 500;                   // Recent window size (number of windows to track)
const double THETA = 200.0;          // Persistence threshold (minimum active windows)
const int DELTA = 100;               // Query step size (query interval in windows)

// ===== Experiment Parameters =====
const int TOTAL_WINDOWS = 0;      // Total number of windows to process (0 = process all)

#endif // PARA_H
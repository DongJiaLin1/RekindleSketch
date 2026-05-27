#ifndef PARM_H
#define PARM_H

#include <cstddef>

//==============================================================================
// Memory Configuration
//==============================================================================

const size_t TOTAL_MEMORY_BYTES = 200 * 1024;  // Total memory budget (200KB)
const int DEFAULT_N = 55;                       // Number of cells per bucket
const size_t ESTIMATED_CELL_SIZE = 41;           // Memory per cell (bytes)

// Calculated number of buckets based on memory constraints
const int CALCULATED_M = static_cast<int>(
    TOTAL_MEMORY_BYTES / (ESTIMATED_CELL_SIZE * DEFAULT_N)
);
const int DEFAULT_M = (CALCULATED_M > 0) ? CALCULATED_M : 1;

//==============================================================================
// Score Decay Parameters
//==============================================================================

const double DEFAULT_ALPHA = 0.015;  // Exponential decay rate
const double DEFAULT_BETA = 1.5;     // Reward function parameter
const double DEFAULT_Q = 1.0;        // Base reward value
const int DEFAULT_Z = 1;            // Replacement probability parameter

//==============================================================================
// Sliding Window Parameters
//==============================================================================

const int R = 500;                   // Number of recent windows to track
const double DALTA = 300.0;          // Persistence threshold (min active windows)
const int S = 100;                   // Query interval (windows between queries)

//==============================================================================
// Experiment Parameters
//==============================================================================

const int TOTAL_WINDOWS = 0;         // Windows to process (0 = process all)
const std::string DATA_FOLDER_PATH = "";  // Path to CSV data files

#endif  // PARM_H

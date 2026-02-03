#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include "RekindleSketch.h"
#include "parm.h"

namespace fs = std::filesystem;

// ===== Data Path Configuration =====
const std::string DATA_FOLDER_PATH = "";  // Path to CSV data files directory

/**
 * Read flow fingerprints from CSV file
 * Expects CSV format with at least 2 columns: IP tuple and fingerprint
 * Skips header row and empty lines
 * @param filepath Path to CSV file
 * @return Vector of flow fingerprint strings
 */
std::vector<std::string> read_csv_flows(const std::string& filepath) {
    std::vector<std::string> flows;
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return flows;
    }
    
    std::string line;
    bool is_first_line = true;  // Skip header row
    
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        
        // Skip header row
        if (is_first_line) {
            is_first_line = false;
            continue;
        }
        
        // Parse CSV line: extract IP tuple (col 1) and fingerprint (col 2)
        std::istringstream iss(line);
        std::string column;
        std::vector<std::string> columns;
        
        // Read all columns
        while (std::getline(iss, column, ',')) {
            // Trim whitespace
            column.erase(0, column.find_first_not_of(" \t"));
            column.erase(column.find_last_not_of(" \t") + 1);
            
            // Remove quotes if present
            if (!column.empty() && column.front() == '"' && column.back() == '"') {
                column = column.substr(1, column.length() - 2);
            }
            
            columns.push_back(column);
        }
        
        // Require at least 2 columns: IP tuple and fingerprint
        if (columns.size() >= 2) {
            std::string fingerprint = columns[1];  // Second column: 16-bit fingerprint
            
            // Trim fingerprint whitespace
            fingerprint.erase(0, fingerprint.find_first_not_of(" \t"));
            fingerprint.erase(fingerprint.find_last_not_of(" \t") + 1);
            
            if (!fingerprint.empty()) {
                flows.push_back(fingerprint);
            }
        }
    }
    
    // File automatically closed by RAII
    return flows;
}

/**
 * Get all CSV file paths from a directory (sorted by filename)
 * @param folder_path Directory path to search
 * @return Vector of CSV file paths, sorted alphabetically
 */
std::vector<std::string> get_csv_files(const std::string& folder_path) {
    std::vector<std::string> csv_files;
    
    try {
        if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
            std::cerr << "Error: Folder does not exist or is not a valid directory: " << folder_path << std::endl;
            return csv_files;
        }
        
        // Iterate through all files in directory
        for (const auto& entry : fs::directory_iterator(folder_path)) {
            if (entry.is_regular_file()) {
                std::string filepath = entry.path().string();
                std::string extension = entry.path().extension().string();
                
                // Convert to lowercase for comparison
                std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
                
                if (extension == ".csv") {
                    csv_files.push_back(filepath);
                }
            }
        }
        
        // Sort files by filename
        std::sort(csv_files.begin(), csv_files.end());
        
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error: Problem accessing folder: " << e.what() << std::endl;
    }
    
    return csv_files;
}

/**
 * Main function: processes CSV files and performs sliding window queries
 * Each CSV file represents one time window
 * Performs periodic queries at DELTA intervals after R windows have been processed
 */
int main() {
    RekindleSketch sketch;
    
    // Get all CSV files from data directory
    std::vector<std::string> csv_files = get_csv_files(DATA_FOLDER_PATH);
    if (csv_files.empty()) {
        std::cerr << "Error: No CSV files found" << std::endl;
        return 1;
    }
    
    int processed_files = 0;
    
    // Process each CSV file as a time window
    for (size_t i = 0; i < csv_files.size() && (TOTAL_WINDOWS == 0 || processed_files < TOTAL_WINDOWS); i++) {
        std::vector<std::string> flows = read_csv_flows(csv_files[i]);
        
        if (flows.empty()) continue;
        
        // Update sketch with all flows in current window
        for (const auto& fingerprint : flows) {
            sketch.update(fingerprint);
        }
        
        processed_files++;
        sketch.end_window();  // End current window and update decay counters
        
        // Perform sliding window query (only executes at DELTA intervals)
        sketch.sliding_window_query(THETA);
    }
    
    // Final query check: ensure last window is queried if needed
    // This handles edge case where loop ends exactly at a query interval
    if (processed_files == TOTAL_WINDOWS && TOTAL_WINDOWS > 0) {
        int current_window = sketch.get_current_window_id();
        int actual_end_window = current_window - 1;
        int query_interval = (actual_end_window - R) % DELTA;
        
        if (current_window > R && query_interval == 0) {
            const auto& sliding_results = sketch.get_sliding_window_results();
            bool already_queried = false;
            if (!sliding_results.empty()) {
                const auto& latest_result = sliding_results.back();
                if (latest_result.end_window == actual_end_window) {
                    already_queried = true;
                }
            }
            
            if (!already_queried) {
                sketch.sliding_window_query(THETA);
            }
        }
    }
    
    return 0;
}

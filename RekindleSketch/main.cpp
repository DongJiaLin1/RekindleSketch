#include <algorithm>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "RekindleSketch.h"
#include "parm.h"

namespace fs = std::filesystem;

//==============================================================================
// Helper Functions
//==============================================================================

static inline std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

static std::vector<std::string> read_csv_flows(const std::string& filepath) {
    std::vector<std::string> flows;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return flows;
    }

    std::string line;
    bool is_header = true;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (is_header) {
            is_header = false;
            continue;
        }

        std::istringstream iss(line);
        std::string column;
        std::vector<std::string> columns;

        while (std::getline(iss, column, ',')) {
            column = trim(column);
            if (!column.empty() && column.front() == '"' && column.back() == '"') {
                column = column.substr(1, column.length() - 2);
            }
            columns.push_back(column);
        }

        if (columns.size() >= 2) {
            std::string fingerprint = trim(columns[1]);
            if (!fingerprint.empty()) {
                flows.push_back(fingerprint);
            }
        }
    }

    return flows;
}

static std::vector<std::string> get_csv_files(const std::string& folder_path) {
    std::vector<std::string> csv_files;

    try {
        if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
            std::cerr << "Error: Invalid folder path: " << folder_path << std::endl;
            return csv_files;
        }

        for (const auto& entry : fs::directory_iterator(folder_path)) {
            if (entry.is_regular_file()) {
                std::string ext = entry.path().extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                if (ext == ".csv") {
                    csv_files.push_back(entry.path().string());
                }
            }
        }

        std::sort(csv_files.begin(), csv_files.end());
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return csv_files;
}

//==============================================================================
// Main Entry Point
//==============================================================================

int main() {
    RekindleSketch sketch;
    std::vector<std::string> csv_files = get_csv_files(DATA_FOLDER_PATH);

    if (csv_files.empty()) {
        std::cerr << "Error: No CSV files found" << std::endl;
        return 1;
    }

    int processed_files = 0;

    for (size_t i = 0; i < csv_files.size() &&
         (TOTAL_WINDOWS == 0 || processed_files < TOTAL_WINDOWS); i++) {
        std::vector<std::string> flows = read_csv_flows(csv_files[i]);
        if (flows.empty()) continue;

        for (const auto& fingerprint : flows) {
            sketch.update(fingerprint);
        }

        processed_files++;
        sketch.end_window();
        sketch.sliding_window_query(DALTA);
    }

    // Handle edge case: query at the end if needed
    if (processed_files == TOTAL_WINDOWS && TOTAL_WINDOWS > 0) {
        int current_window = sketch.get_current_window_id();
        int end_window = current_window - 1;
        int query_interval = (end_window - R) % S;

        if (current_window > R && query_interval == 0) {
            const auto& results = sketch.get_sliding_window_results();
            bool already_queried = !results.empty() &&
                                   results.back().end_window == end_window;
            if (!already_queried) {
                sketch.sliding_window_query(DALTA);
            }
        }
    }

    return 0;
}

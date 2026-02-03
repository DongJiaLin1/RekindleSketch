#ifndef REKINDLESKETCH_H
#define REKINDLESKETCH_H

#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <memory>
#include <deque>
#include "parm.h"
#include "MurmurHash3.h"

/**
 * RekindleSketch: A memory-efficient data structure for tracking persistent flows
 * Implements a sketch-based algorithm with exponential decay and reward mechanisms
 */
class RekindleSketch {
private:
    /**
     * Cell: Represents a single memory cell in the sketch
     * Stores flow fingerprint, activity flag, memory score, and decay information
     */
    struct Cell {
        uint16_t key;                    // Flow fingerprint (0 = empty)
        bool flag;                       // Current window active flag
        double score;                    // Memory score (decays exponentially)
        uint16_t decay;                 // Missing windows count (decay factor)
        uint16_t last_active_window;     // Last active window ID

        Cell() : key(0), flag(false), score(0.0), decay(0), last_active_window(-1) {}
        bool is_empty() const { return key == 0; }
        void clear() {
            key = 0;
            flag = false;
            score = 0.0;
            decay = 0;
            last_active_window = -1;
        }
    };

    /**
     * Bucket: Container for multiple cells
     * Provides cell lookup, replacement candidate selection, and management operations
     */
    struct Bucket {
        int bucket_id;
        std::vector<Cell> cells;

        Bucket(int id, int n) : bucket_id(id), cells(n) {}

        Cell& get_cell(int k) {
            if (k < 1 || k > static_cast<int>(cells.size())) {
                throw std::out_of_range("Cell index out of range");
            }
            return cells[k-1];
        }

        const Cell& get_cell(int k) const {
            if (k < 1 || k > static_cast<int>(cells.size())) {
                throw std::out_of_range("Cell index out of range");
            }
            return cells[k-1];
        }

        int find_cell_index(uint16_t flow_id) const {
            for (size_t i = 0; i < cells.size(); i++) {
                if (!cells[i].is_empty() && cells[i].key == flow_id) {
                    return static_cast<int>(i) + 1;
                }
            }
            return 0;
        }

        int find_empty_cell_index() const {
            for (size_t i = 0; i < cells.size(); i++) {
                if (cells[i].is_empty()) {
                    return static_cast<int>(i) + 1;
                }
            }
            return 0;
        }

        /**
         * Find replacement candidate: selects cell with maximum decay value
         * among inactive cells (cells not active in current window)
         */
        int find_replacement_candidate_index() const {
            int candidate_idx = 0;
            int max_decay = -1;

            for (size_t i = 0; i < cells.size(); i++) {
                const Cell& cell = cells[i];
                if (!cell.is_empty() && !cell.flag) {
                    if (cell.decay > max_decay) {
                        max_decay = cell.decay;
                        candidate_idx = static_cast<int>(i) + 1;
                    }
                }
            }

            return candidate_idx;
        }
    };

    // Sketch configuration parameters
    const int m_buckets;              // Number of hash buckets
    const int m_cells_per_bucket;     // Number of cells per bucket
    const double m_alpha;             // Exponential decay rate
    const double m_beta;              // Reward function parameter
    const double m_Q;                 // Base reward value
    const int m_Z;                    // Replacement probability parameter

    std::vector<std::unique_ptr<Bucket>> buckets;
    uint32_t hash_seed;               // Seed for hash function

    // Window management
    int current_window_id;            // Current window identifier
    int total_windows_processed;      // Total number of windows processed
    
    /**
     * SlidingWindowResult: Stores query results for a sliding window
     * Contains window range and list of persistent flows found
     */
    struct SlidingWindowResult {
        int start_window;
        int end_window;
        std::vector<std::string> persistent_flows;
        
        SlidingWindowResult(int start, int end) 
            : start_window(start), end_window(end) {}
    };
    
    std::deque<SlidingWindowResult> sliding_window_queue;  // Queue for sliding window results

    // Random number generation for probabilistic replacement
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;

    double score_upper_bound;         // Upper bound for memory score: S_max = 2Q / (1 - e^(-α))

public:
    RekindleSketch()
            : m_buckets(DEFAULT_M),
              m_cells_per_bucket(DEFAULT_N),
              m_alpha(DEFAULT_ALPHA),
              m_beta(DEFAULT_BETA),
              m_Q(DEFAULT_Q),
              m_Z(DEFAULT_Z),
              hash_seed(123456),
              current_window_id(1),
              total_windows_processed(0) {
        initialize();
    }

    RekindleSketch(int custom_m, int custom_n,
                   double custom_alpha = DEFAULT_ALPHA,
                   double custom_beta = DEFAULT_BETA,
                   double custom_Q = DEFAULT_Q,
                   int custom_Z = DEFAULT_Z,
                   uint32_t seed = 123456)
            : m_buckets(custom_m),
              m_cells_per_bucket(custom_n),
              m_alpha(custom_alpha),
              m_beta(custom_beta),
              m_Q(custom_Q),
              m_Z(custom_Z),
              hash_seed(seed),
              current_window_id(1),
              total_windows_processed(0) {
        initialize();
    }

    RekindleSketch(int custom_m, size_t total_memory_bytes, size_t cell_size = 16,
                   double custom_alpha = DEFAULT_ALPHA,
                   double custom_beta = DEFAULT_BETA,
                   double custom_Q = DEFAULT_Q,
                   int custom_Z = DEFAULT_Z,
                   uint32_t seed = 123456)
            : m_buckets(custom_m),
              m_cells_per_bucket((custom_m > 0 && cell_size > 0) 
                                 ? std::max(1, static_cast<int>(total_memory_bytes / (cell_size * custom_m)))
                                 : 1),
              m_alpha(custom_alpha),
              m_beta(custom_beta),
              m_Q(custom_Q),
              m_Z(custom_Z),
              hash_seed(seed),
              current_window_id(1),
              total_windows_processed(0) {
        initialize();
    }

private:
    /**
     * Initialize sketch: set up random number generator, calculate score bounds,
     * and allocate bucket memory
     */
    void initialize() {
        std::random_device rd;
        rng = std::mt19937(rd());
        dist = std::uniform_real_distribution<double>(0.0, 1.0);

        // Score upper bound: S_max = 2Q / (1 - e^(-α))
        score_upper_bound = 2.0 * m_Q / (1.0 - std::exp(-m_alpha));

        buckets.reserve(m_buckets);
        for (int i = 1; i <= m_buckets; i++) {
            buckets.push_back(std::make_unique<Bucket>(i, m_cells_per_bucket));
        }
    }

    uint32_t murmur_hash(const std::string& key) const {
        uint32_t hash_value;
        MurmurHash3_x86_32(key.data(), static_cast<int>(key.length()), hash_seed, &hash_value);
        return hash_value;
    }

    /**
     * Convert string to 16-bit fingerprint
     * Tries to parse as number first, falls back to hash if parsing fails
     * @param key Input string
     * @return 16-bit fingerprint (never returns 0)
     */
    uint16_t string_to_fingerprint(const std::string& key) const {
        try {
            unsigned long value = std::stoul(key);
            uint16_t fingerprint = static_cast<uint16_t>(value & 0xFFFF);
            return (fingerprint == 0) ? 1 : fingerprint;
        } catch (const std::exception&) {
            uint32_t hash_value = murmur_hash(key);
            uint16_t fingerprint = static_cast<uint16_t>(hash_value & 0xFFFF);
            return (fingerprint == 0) ? 1 : fingerprint;
        }
    }

    /**
     * Compute bucket index from fingerprint using multiplicative hash
     * @param fingerprint 16-bit flow fingerprint
     * @return Bucket index (0-based)
     */
    int get_bucket_index(uint16_t fingerprint) const {
        uint32_t hash_value = static_cast<uint32_t>(fingerprint) * 2654435761U;
        return hash_value % m_buckets;
    }

public:
    /**
     * Reward function: G(D) = 1 + exp(-β * D)
     * Returns higher reward for flows that reappear after shorter gaps
     * @param D Number of missing windows
     * @return Reward value
     */
    double reward_function(int D) const {
        return 1.0 + std::exp(-m_beta * D);
    }

    /**
     * Calculate memory score: S_W = S_Ŵ * exp(-α * (D + 1)) + Q * G(D)
     * Combines exponential decay with reward mechanism
     * @param prev_score Previous memory score
     * @param D Number of missing windows
     * @param current_window Current window ID (unused, kept for API consistency)
     * @return Updated memory score (capped at score_upper_bound)
     */
    double calculate_memory_score(double prev_score, int D, int /* current_window */) const {
        double decay_factor = std::exp(-m_alpha * (D + 1));
        double reward = m_Q * reward_function(D);
        double new_score = prev_score * decay_factor + reward;
        return std::min(new_score, score_upper_bound);
    }

    /**
     * Update sketch with a flow: convert string ID to fingerprint and update
     * @param flow_id Flow identifier (string format)
     * @return true if update successful, false otherwise
     */
    bool update(const std::string& flow_id) {
        uint16_t fingerprint = string_to_fingerprint(flow_id);
        if (fingerprint == 0) {
            return false;
        }
        return update_by_fingerprint(fingerprint);
    }
    
    /**
     * Update sketch by fingerprint: implements three-step update strategy
     * 1. If flow exists and inactive: update score and activate
     * 2. If empty cell exists: insert new flow
     * 3. If bucket full: probabilistic replacement based on decay
     * @param fingerprint 16-bit flow fingerprint
     * @return true if update successful, false otherwise
     */
    bool update_by_fingerprint(uint16_t fingerprint) {
        int bucket_idx = get_bucket_index(fingerprint);
        Bucket& bucket = *buckets[bucket_idx];

        int cell_idx = bucket.find_cell_index(fingerprint);

        if (cell_idx > 0) {
            Cell& cell = bucket.get_cell(cell_idx);
            if (!cell.flag) {
                cell.score = calculate_memory_score(cell.score, cell.decay, current_window_id);
                cell.flag = true;
                cell.decay = 0;
                cell.last_active_window = current_window_id;
            }
            return true;
        }

        int empty_cell_idx = bucket.find_empty_cell_index();
        if (empty_cell_idx > 0) {
            Cell& cell = bucket.get_cell(empty_cell_idx);
            cell.key = fingerprint;
            cell.flag = true;
            cell.score = m_Q;
            cell.decay = 0;
            cell.last_active_window = current_window_id;
            return true;
        }

        // Step 3: Probabilistic replacement if bucket is full
        int candidate_idx = bucket.find_replacement_candidate_index();
        if (candidate_idx > 0) {
            Cell& candidate = bucket.get_cell(candidate_idx);
            double decayed_score = candidate.score * std::exp(-m_alpha * candidate.decay);
            // Replacement probability: P = Z / (Z + decayed_score)
            double replace_prob = m_Z / (m_Z + decayed_score);

            if (dist(rng) < replace_prob) {
                candidate.clear();
                candidate.key = fingerprint;
                candidate.flag = true;
                candidate.score = m_Q;
                candidate.decay = 0;
                candidate.last_active_window = current_window_id;
                return true;
            }
            return false;
        }

        return false;
    }

    /**
     * Update sketch with all flows in current window
     * @param flows Vector of flow identifiers
     */
    void update_window(const std::vector<std::string>& flows) {
        for (const auto& flow : flows) {
            update(flow);
        }
    }

    /**
     * End current window: update decay counters and reset activity flags
     * Increments decay for inactive cells and advances window counter
     */
    void end_window() {
        for (int d = 1; d <= m_buckets; d++) {
            Bucket& bucket = *buckets[d-1];
            for (int k = 1; k <= m_cells_per_bucket; k++) {
                Cell& cell = bucket.get_cell(k);
                if (!cell.is_empty() && !cell.flag) {
                    cell.decay++;
                }
                cell.flag = false;
            }
        }
        current_window_id++;
        total_windows_processed++;
    }

    /**
     * Query memory score: returns decayed memory score for a flow
     * Calculates latest score considering exponential decay: S = score * exp(-α * decay)
     * @param flow_id Flow identifier (string format)
     * @return Latest memory score, or 0 if flow not found
     */
    double query(const std::string& flow_id) {
        uint16_t fingerprint = string_to_fingerprint(flow_id);
        if (fingerprint == 0) {
            return 0.0;  // Invalid fingerprint
        }
        
        return query_by_fingerprint(fingerprint);
    }

    /**
     * Query by fingerprint: internal query implementation
     * @param fingerprint 16-bit flow fingerprint
     * @return Decayed memory score, or 0 if flow not found
     */
    double query_by_fingerprint(uint16_t fingerprint) {
        int bucket_idx = get_bucket_index(fingerprint);
        const Bucket& bucket = *buckets[bucket_idx];

        int cell_idx = bucket.find_cell_index(fingerprint);

        if (cell_idx > 0) {
            const Cell& cell = bucket.get_cell(cell_idx);
            double latest_score = cell.score * std::exp(-m_alpha * cell.decay);
            return latest_score;
        }

        return 0.0;
    }

    /**
     * Query recent persistence: estimates number of active windows in recent R windows
     * Normalizes memory score and scales by R to estimate persistence
     * @param flow_id Flow identifier (string format)
     * @return Estimated number of active windows, or 0 if flow not found
     */
    double query_recent_persistence(const std::string& flow_id) const {
        uint16_t fingerprint = string_to_fingerprint(flow_id);
        if (fingerprint == 0) {
            return 0.0;
        }
        
        return query_recent_persistence_by_fingerprint(fingerprint);
    }
    
    /**
     * Query recent persistence by fingerprint: internal implementation
     * @param fingerprint 16-bit flow fingerprint
     * @return Estimated active windows = (decayed_score / S_max) * R
     */
    double query_recent_persistence_by_fingerprint(uint16_t fingerprint) const {
        int bucket_idx = get_bucket_index(fingerprint);
        const Bucket& bucket = *buckets[bucket_idx];

        int cell_idx = bucket.find_cell_index(fingerprint);

        if (cell_idx > 0) {
            const Cell& cell = bucket.get_cell(cell_idx);

            double recent_score = cell.score * std::exp(-m_alpha * cell.decay);

            double normalized_score = recent_score / score_upper_bound;

            double estimated_active_windows = normalized_score * R;

            return estimated_active_windows;
        }

        return 0.0;
    }

    /**
     * Find recent persistent flows: scans all cells and returns flows with
     * persistence above threshold
     * @param threshold Persistence threshold (default: THETA)
     * @return Vector of persistent flow identifiers
     */
    std::vector<std::string> find_recent_persistent_flows(double threshold = THETA) {
        std::vector<std::string> persistent_flows;

        for (int d = 1; d <= m_buckets; d++) {
            const Bucket& bucket = *buckets[d-1];

            for (int k = 1; k <= m_cells_per_bucket; k++) {
                const Cell& cell = bucket.get_cell(k);
                if (!cell.is_empty()) {
                    double persistence = query_recent_persistence_by_fingerprint(cell.key);
                    if (persistence >= threshold) {
                        persistent_flows.push_back(std::to_string(cell.key));
                    }
                }
            }
        }

        return persistent_flows;
    }

    /**
     * Sliding window query: performs periodic query every DELTA windows
     * Queries are only performed when current_window_id >= R and at DELTA intervals
     * Results are stored in sliding_window_queue with max size R/DELTA
     * @param threshold Persistence threshold (default: THETA)
     * @return true if query was performed, false otherwise
     */
    bool sliding_window_query(double threshold = THETA) {
        if (current_window_id < R) {
            return false;
        }

        int actual_end_window = current_window_id - 1;
        int start_window = actual_end_window - R + 1;
        int end_window = actual_end_window;
        int query_interval = (actual_end_window - R) % DELTA;
        if (query_interval != 0) {
            return false;
        }

        std::vector<std::string> persistent_flows = find_recent_persistent_flows(threshold);

        SlidingWindowResult result(start_window, end_window);
        result.persistent_flows = persistent_flows;

        sliding_window_queue.push_back(result);

        int max_queue_size = R / DELTA;
        if (sliding_window_queue.size() > max_queue_size) {
            sliding_window_queue.pop_front();
        }
        
        return true;
    }

    const std::deque<SlidingWindowResult>& get_sliding_window_results() const {
        return sliding_window_queue;
    }

    void clear_sliding_window_queue() {
        sliding_window_queue.clear();
    }

    int get_current_window_id() const {
        return current_window_id;
    }

    int get_total_windows_processed() const {
        return total_windows_processed;
    }

    double get_score_upper_bound() const {
        return score_upper_bound;
    }

    /**
     * Reset sketch: clears all cells and resets window counters
     * Useful for reusing sketch instance with new data stream
     */
    void reset() {
        current_window_id = 1;
        total_windows_processed = 0;

        // Clear all cells
        for (int d = 1; d <= m_buckets; d++) {
            Bucket& bucket = *buckets[d-1];

            for (int k = 1; k <= m_cells_per_bucket; k++) {
                Cell& cell = bucket.get_cell(k);
                cell.clear();
            }
        }
    }

    const Bucket& get_bucket(int d) const {
        if (d < 1 || d > m_buckets) {
            throw std::out_of_range("Bucket index out of range");
        }
        return *buckets[d-1];
    }

    const Cell& get_cell(int d, int k) const {
        const Bucket& bucket = get_bucket(d);
        return bucket.get_cell(k);
    }

};

#endif // REKINDLESKETCH_H

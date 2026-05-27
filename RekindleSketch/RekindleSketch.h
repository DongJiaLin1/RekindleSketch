#ifndef REKINDLESKETCH_H
#define REKINDLESKETCH_H

#include <algorithm>
#include <cmath>
#include <deque>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "MurmurHash3.h"
#include "parm.h"



class RekindleSketch {
public:
    // Constructor with memory-based configuration
    RekindleSketch(int custom_m, size_t total_memory_bytes, size_t cell_size = 16,
                   double custom_alpha = DEFAULT_ALPHA,
                   double custom_beta = DEFAULT_BETA,
                   double custom_Q = DEFAULT_Q,
                   int custom_Z = DEFAULT_Z,
                   uint32_t seed = 123456);

    // Constructor with explicit bucket/cell configuration
    RekindleSketch(int custom_m, int custom_n,
                   double custom_alpha = DEFAULT_ALPHA,
                   double custom_beta = DEFAULT_BETA,
                   double custom_Q = DEFAULT_Q,
                   int custom_Z = DEFAULT_Z,
                   uint32_t seed = 123456);

    // Default constructor
    RekindleSketch();

    // Update sketch with a flow
    bool update(const std::string& flow_id);
    bool update_by_fingerprint(uint16_t fingerprint);

    // Update multiple flows for current window
    void update_window(const std::vector<std::string>& flows);

    // End current window and apply decay
    void end_window();

    // Query functions
    double query(const std::string& flow_id);
    double query_by_fingerprint(uint16_t fingerprint);
    double query_persistence(const std::string& flow_id) const;
    double query_persistence_by_fingerprint(uint16_t fingerprint) const;

    // Find all persistent flows above threshold
    std::vector<std::string> find_persistent_flows(double threshold = DALTA);

    // Sliding window query (periodic)
    bool sliding_window_query(double threshold = DALTA);
    const std::deque<struct SlidingWindowResult>& get_sliding_window_results() const;

    // State accessors
    int current_window() const { return current_window_id_; }
    int total_windows_processed() const { return total_windows_processed_; }

    // Reset sketch state
    void reset();

private:

    struct Cell {
        uint16_t key = 0;
        bool flag = false;
        uint8_t score = 0;
        uint8_t decay = 0;
        int last_active_window = -1;

        bool empty() const { return key == 0; }
        void clear() { key = 0; flag = false; score = 0; decay = 0; last_active_window = -1; }
    };

    struct Bucket {
        int id;
        std::vector<Cell> cells;

        explicit Bucket(int bucket_id, int num_cells) : id(bucket_id), cells(num_cells) {}

        Cell& at(size_t idx) { return cells.at(idx); }
        const Cell& at(size_t idx) const { return cells.at(idx); }

        int find_cell(uint16_t flow_id) const {
            for (size_t i = 0; i < cells.size(); ++i) {
                if (!cells[i].empty() && cells[i].key == flow_id) {
                    return static_cast<int>(i);
                }
            }
            return -1;
        }

        int find_empty() const {
            for (size_t i = 0; i < cells.size(); ++i) {
                if (cells[i].empty()) return static_cast<int>(i);
            }
            return -1;
        }

        int find_replacement_candidate() const {
            int candidate = -1;
            int max_decay = -1;
            for (size_t i = 0; i < cells.size(); ++i) {
                if (!cells[i].empty() && !cells[i].flag) {
                    if (static_cast<int>(cells[i].decay) > max_decay) {
                        max_decay = cells[i].decay;
                        candidate = static_cast<int>(i);
                    }
                }
            }
            return candidate;
        }
    };

    struct SlidingWindowResult {
        int start_window;
        int end_window;
        std::vector<std::string> persistent_flows;

        SlidingWindowResult(int start, int end) : start_window(start), end_window(end) {}
    };


    const int num_buckets_;
    const int cells_per_bucket_;
    const double alpha_;
    const double beta_;
    const double Q_;
    const int Z_;

    std::vector<std::unique_ptr<Bucket>> buckets_;
    uint32_t hash_seed_;
    int current_window_id_ = 1;
    int total_windows_processed_ = 0;

    std::deque<SlidingWindowResult> sliding_window_results_;
    std::mt19937 rng_;
    std::uniform_real_distribution<double> dist_;
    double score_upper_bound_;


    void initialize();

    uint32_t murmur_hash(const std::string& key) const {
        uint32_t hash_value;
        MurmurHash3_x86_32(key.data(), static_cast<int>(key.length()), hash_seed_, &hash_value);
        return hash_value;
    }

    uint16_t string_to_fingerprint(const std::string& key) const {
        try {
            unsigned long value = std::stoul(key);
            uint16_t fp = static_cast<uint16_t>(value & 0xFFFF);
            return (fp == 0) ? 1 : fp;
        } catch (const std::exception&) {
            uint32_t hash_value = murmur_hash(key);
            uint16_t fp = static_cast<uint16_t>(hash_value & 0xFFFF);
            return (fp == 0) ? 1 : fp;
        }
    }

    int get_bucket_index(uint16_t fingerprint) const {
        uint32_t hash_value = static_cast<uint32_t>(fingerprint) * 2654435761U;
        return static_cast<int>(hash_value % num_buckets_);
    }

    double reward_function(int decay) const {
        return 1.0 + std::exp(-beta_ * decay);
    }

    uint8_t calculate_score(uint8_t prev_score, uint8_t decay) const {
        double prev = static_cast<double>(prev_score);
        double decay_factor = std::exp(-alpha_ * (decay + 1));
        double new_score = prev * decay_factor + Q_ * reward_function(decay);
        return static_cast<uint8_t>(std::min(new_score, 255.0));
    }
};


inline RekindleSketch::RekindleSketch()
    : RekindleSketch(DEFAULT_M, DEFAULT_N) {}

inline RekindleSketch::RekindleSketch(int custom_m, int custom_n,
                                       double custom_alpha, double custom_beta,
                                       double custom_Q, int custom_Z, uint32_t seed)
    : num_buckets_(custom_m),
      cells_per_bucket_(custom_n),
      alpha_(custom_alpha),
      beta_(custom_beta),
      Q_(custom_Q),
      Z_(custom_Z),
      hash_seed_(seed) {
    initialize();
}

inline RekindleSketch::RekindleSketch(int custom_m, size_t total_memory_bytes,
                                       size_t cell_size, double custom_alpha,
                                       double custom_beta, double custom_Q,
                                       int custom_Z, uint32_t seed)
    : num_buckets_(custom_m),
      cells_per_bucket_(custom_m > 0 && cell_size > 0
                         ? std::max(1, static_cast<int>(total_memory_bytes / (cell_size * custom_m)))
                         : 1),
      alpha_(custom_alpha),
      beta_(custom_beta),
      Q_(custom_Q),
      Z_(custom_Z),
      hash_seed_(seed) {
    initialize();
}

inline void RekindleSketch::initialize() {
    std::random_device rd;
    rng_ = std::mt19937(rd());
    dist_ = std::uniform_real_distribution<double>(0.0, 1.0);

    score_upper_bound_ = 2.0 * Q_ / (1.0 - std::exp(-alpha_));

    buckets_.reserve(num_buckets_);
    for (int i = 0; i < num_buckets_; ++i) {
        buckets_.push_back(std::make_unique<Bucket>(i, cells_per_bucket_));
    }
}

inline bool RekindleSketch::update(const std::string& flow_id) {
    uint16_t fingerprint = string_to_fingerprint(flow_id);
    if (fingerprint == 0) return false;
    return update_by_fingerprint(fingerprint);
}

inline bool RekindleSketch::update_by_fingerprint(uint16_t fingerprint) {
    int bucket_idx = get_bucket_index(fingerprint);
    Bucket& bucket = *buckets_[bucket_idx];

    int cell_idx = bucket.find_cell(fingerprint);
    if (cell_idx >= 0) {
        Cell& cell = bucket.at(static_cast<size_t>(cell_idx));
        if (!cell.flag) {
            cell.score = calculate_score(cell.score, cell.decay);
            cell.flag = true;
            cell.decay = 0;
            cell.last_active_window = current_window_id_;
        }
        return true;
    }

    cell_idx = bucket.find_empty();
    if (cell_idx >= 0) {
        Cell& cell = bucket.at(static_cast<size_t>(cell_idx));
        cell.key = fingerprint;
        cell.flag = true;
        cell.score = static_cast<uint8_t>(Q_);
        cell.decay = 0;
        cell.last_active_window = current_window_id_;
        return true;
    }

    // Probabilistic replacement
    cell_idx = bucket.find_replacement_candidate();
    if (cell_idx >= 0) {
        Cell& candidate = bucket.at(static_cast<size_t>(cell_idx));
        double decayed_score = static_cast<double>(candidate.score) *
                               std::exp(-alpha_ * candidate.decay);
        double replace_prob = static_cast<double>(Z_) /
                             (Z_ + decayed_score);

        if (dist_(rng_) < replace_prob) {
            candidate.clear();
            candidate.key = fingerprint;
            candidate.flag = true;
            candidate.score = static_cast<uint8_t>(Q_);
            candidate.decay = 0;
            candidate.last_active_window = current_window_id_;
            return true;
        }
    }

    return false;
}

inline void RekindleSketch::update_window(const std::vector<std::string>& flows) {
    for (const auto& flow : flows) {
        update(flow);
    }
}

inline void RekindleSketch::end_window() {
    for (auto& bucket_ptr : buckets_) {
        for (auto& cell : bucket_ptr->cells) {
            if (!cell.empty() && !cell.flag) {
                if (cell.decay < 255) cell.decay++;
            }
            cell.flag = false;
        }
    }
    ++current_window_id_;
    ++total_windows_processed_;
}

inline double RekindleSketch::query(const std::string& flow_id) {
    uint16_t fingerprint = string_to_fingerprint(flow_id);
    if (fingerprint == 0) return 0.0;
    return query_by_fingerprint(fingerprint);
}

inline double RekindleSketch::query_by_fingerprint(uint16_t fingerprint) {
    int bucket_idx = get_bucket_index(fingerprint);
    const Bucket& bucket = *buckets_[bucket_idx];

    int cell_idx = bucket.find_cell(fingerprint);
    if (cell_idx >= 0) {
        const Cell& cell = bucket.at(static_cast<size_t>(cell_idx));
        return static_cast<double>(cell.score) *
               std::exp(-alpha_ * cell.decay);
    }
    return 0.0;
}

inline double RekindleSketch::query_persistence(const std::string& flow_id) const {
    uint16_t fingerprint = string_to_fingerprint(flow_id);
    if (fingerprint == 0) return 0.0;
    return query_persistence_by_fingerprint(fingerprint);
}

inline double RekindleSketch::query_persistence_by_fingerprint(uint16_t fingerprint) const {
    int bucket_idx = get_bucket_index(fingerprint);
    const Bucket& bucket = *buckets_[bucket_idx];

    int cell_idx = bucket.find_cell(fingerprint);
    if (cell_idx >= 0) {
        const Cell& cell = bucket.at(static_cast<size_t>(cell_idx));
        double recent_score = static_cast<double>(cell.score) *
                             std::exp(-alpha_ * cell.decay);
        double normalized_score = recent_score / score_upper_bound_;
        return normalized_score * R;
    }
    return 0.0;
}

inline std::vector<std::string> RekindleSketch::find_persistent_flows(double threshold) {
    std::vector<std::string> persistent_flows;

    for (const auto& bucket_ptr : buckets_) {
        for (const auto& cell : bucket_ptr->cells) {
            if (!cell.empty()) {
                double persistence = query_persistence_by_fingerprint(cell.key);
                if (persistence >= threshold) {
                    persistent_flows.push_back(std::to_string(cell.key));
                }
            }
        }
    }

    return persistent_flows;
}

inline bool RekindleSketch::sliding_window_query(double threshold) {
    if (current_window_id_ <= R) return false;

    int end_window = current_window_id_ - 1;
    int start_window = end_window - R + 1;
    int query_interval = (end_window - R) % S;
    if (query_interval != 0) return false;

    std::vector<std::string> flows = find_persistent_flows(threshold);
    sliding_window_results_.push_back(SlidingWindowResult(start_window, end_window));
    sliding_window_results_.back().persistent_flows = std::move(flows);

    size_t max_size = static_cast<size_t>(R / S);
    while (sliding_window_results_.size() > max_size) {
        sliding_window_results_.pop_front();
    }

    return true;
}

inline const std::deque<struct RekindleSketch::SlidingWindowResult>&
RekindleSketch::get_sliding_window_results() const {
    return sliding_window_results_;
}

inline void RekindleSketch::reset() {
    current_window_id_ = 1;
    total_windows_processed_ = 0;
    sliding_window_results_.clear();

    for (auto& bucket_ptr : buckets_) {
        for (auto& cell : bucket_ptr->cells) {
            cell.clear();
        }
    }
}

#endif  // REKINDLESKETCH_H

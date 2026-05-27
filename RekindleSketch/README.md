# RekindleSketch

A memory-efficient sketch data structure for detecting persistent network flows in sliding windows.

## Overview

RekindleSketch is a time-aware sketch-based data structure for detecting recent persistent flows that remain continuously active within the latest R windows. It introduces a memory score mechanism to enable recent persistent flows to be progressively reinforced while obsolete historical activity is quickly forgotten. Based on this mechanism, RekindleSketch further employs a probabilistic replacement strategy to preserve recent persistent flows under hash collisions and efficiently evict inactive flows, achieving high accuracy with compact memory usage.

## Features

- **Memory-efficient**: Configurable memory constraints with automatic bucket/cell calculation
- **Time-aware scoring:**: Memory-score-based tracking with decay and reward mechanisms for recent persistent flows
- **Configurable parameters**: Adjustable decay coefficients, reward functions, and thresholds
- **Fast**: Optimized data structures and inline implementations

## Requirements

- C++17 or later
- CMake 3.15+
- Compiler: GCC 7+, Clang 5+, MSVC 2017+

## Build

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

### Windows (Visual Studio)

```bash
mkdir build && cd build
cmake .. -G "Visual Studio 17 2022"
cmake --build . --config Release
```

## Configuration

Edit `parm.h` to configure:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TOTAL_MEMORY_BYTES` | 200KB | Total memory budget |
| `DEFAULT_M` | auto | Number of hash buckets |
| `DEFAULT_N` | 55 | Cells per bucket |
| `DEFAULT_ALPHA` | 0.015 | Exponential decay factor |
| `DEFAULT_BETA` | 1.5 | Reward coefficient |
| `R` | 500 | Window size |
| `DELTA` | 300 | Persistence threshold |


## Usage

1. Set `DATA_FOLDER_PATH` in `parm.h` to your CSV data directory
2. CSV format: first column = IP tuple, second column = fingerprint
3. Build and run:

```bash
./RekindleSketch
```

## Project Structure

```
.
├── CMakeLists.txt
├── main.cpp
├── RekindleSketch.h
├── parm.h
├── MurmurHash3.h
├── README.md
└── LICENSE
```

## Algorithm

RekindleSketch uses:
- **Memory score**: Tracks flow activity with exponential decay
- **Time-awareness**: Detects recent persistent flows over the latest R windows
- **Probabilistic replacement**: Handles hash collisions with p = \frac{1}{\max(1, Z - R_{h(f)k}.D) \times \hat{S} + 1}

## License

MIT License

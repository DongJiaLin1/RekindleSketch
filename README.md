# RekindleSketch

A memory-efficient sketch data structure for detecting persistent network flows in sliding windows.

## Overview

RekindleSketch implements a sketch-based algorithm for identifying persistent flows that remain active across multiple time windows. It uses exponential decay for score tracking and probabilistic replacement for handling hash collisions.

## Features

- **Memory-efficient**: Configurable memory constraints with automatic bucket/cell calculation
- **Sliding window queries**: Efficient querying of persistent flows in recent windows
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
| `DEFAULT_N` | 550 | Cells per bucket |
| `DEFAULT_ALPHA` | 0.015 | Decay rate |
| `R` | 500 | Sliding window size |
| `THETA` | 200.0 | Persistence threshold |
| `DELTA` | 100 | Query interval |

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
- **Sliding window**: Queries persistent flows in recent R windows
- **Probabilistic replacement**: Handles hash collisions with P = Z / (Z + S_decayed)

## License

MIT License

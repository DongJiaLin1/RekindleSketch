RekindleSketch
RekindleSketch is a memory-efficient data structure for detecting and tracking persistent network flows in sliding windows.

Overview
RekindleSketch implements a sketch-based algorithm for identifying persistent flows that remain active across multiple time windows. It uses a memory-efficient design with configurable memory constraints and provides accurate persistence estimation.

Features
Memory-efficient: Configurable memory usage with automatic bucket/cell calculation
Sliding window queries: Efficient querying of persistent flows in recent windows
Configurable parameters: Adjustable decay coefficients, reward functions, and thresholds
High performance: Optimized data structures and algorithms
Requirements
C++17 or later
CMake 3.15 or later
Compiler supporting C++17 (GCC 7+, Clang 5+, MSVC 2017+)
Building
Linux/macOS
mkdir build
cd build
cmake ..
cmake --build .
Windows (Visual Studio)
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022"
cmake --build . --config Release
Configuration
Edit parm.h to configure:

Usage
Configure DATA_FOLDER_PATH in parm.h to point to your CSV files directory
Ensure CSV files follow the required format (see Data Format below)
Build and run:
./RekindleSketch_1_14
Data Format
The program expects CSV files with:

First column: IP address tuple (optional, not used in current implementation)
Second column: 16-bit fingerprint (required)

Project Structure
.
├── CMakeLists.txt      # CMake build configuration
├── main.cpp            # Main program entry point
├── RekindleSketch.h    # Core RekindleSketch implementation
├── parm.h              # Configuration parameters
├── MurmurHash3.h       # Hash function implementation
├── README.md           # This file
├── LICENSE             # License information
└── .gitignore          # Git ignore rules
Algorithm
RekindleSketch uses:

Memory score: Tracks flow activity with exponential decay
Sliding window: Queries persistent flows in recent R windows
Probabilistic replacement: Handles hash collisions efficiently
License
MIT License - see LICENSE file for details.

Contributing
Contributions are welcome! Please ensure code follows existing style and includes appropriate tests.

References
This implementation is based on the RekindleSketch algorithm for persistent flow detection in network monitoring.

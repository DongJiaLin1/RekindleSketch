cmake_minimum_required(VERSION 3.15)
project(RekindleSketch VERSION 1.14 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Source files
set(SOURCES
    main.cpp
    RekindleSketch.h
    parm.h
    MurmurHash3.h
)

add_executable(RekindleSketch ${SOURCES})

# Link filesystem library (required for C++17 on GCC/Clang)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_link_libraries(RekindleSketch PUBLIC stdc++fs)
endif()

# Set output directory
set_target_properties(RekindleSketch PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

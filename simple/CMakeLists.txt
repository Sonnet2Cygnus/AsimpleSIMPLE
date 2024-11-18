cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(simpleSolver
        VERSION 1.0
        LANGUAGES CXX
        )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add paths for head files
include_directories(.)

link_directories()

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(simpleSolver ${SRC})

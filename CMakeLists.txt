cmake_minimum_required(VERSION 3.14)
project(BioMesh2 VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Include directories
include_directories(include)

# Create library for the core functionality
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_library(biomesh2_lib ${SOURCES})
target_include_directories(biomesh2_lib PUBLIC include)

# Create executable for example usage
add_executable(biomesh2_example examples/main.cpp)
target_link_libraries(biomesh2_example biomesh2_lib)

# Create executable for octree demonstration
add_executable(octree_demo examples/octree_demo.cpp)
target_link_libraries(octree_demo biomesh2_lib)

# Create executable for mesh generation demonstration
add_executable(mesh_demo examples/mesh_demo.cpp)
target_link_libraries(mesh_demo biomesh2_lib)

# Create executable for mesh generation tests
add_executable(mesh_test examples/mesh_test.cpp)
target_link_libraries(mesh_test biomesh2_lib)

# Create executable for comprehensive mesh demo
add_executable(comprehensive_demo examples/comprehensive_demo.cpp)
target_link_libraries(comprehensive_demo biomesh2_lib)

# Create executable for GiD export demonstration
add_executable(gid_export_demo examples/gid_export_demo.cpp)
target_link_libraries(gid_export_demo biomesh2_lib)

# Create executable for complete pipeline demonstration
add_executable(complete_pipeline_demo examples/complete_pipeline_demo.cpp)
target_link_libraries(complete_pipeline_demo biomesh2_lib)

# Enable testing
enable_testing()

# Find GoogleTest
find_package(GTest QUIET)
if(GTest_FOUND)
    # Create test executable
    file(GLOB_RECURSE TEST_SOURCES "tests/*.cpp")
    add_executable(biomesh2_tests ${TEST_SOURCES})
    target_link_libraries(biomesh2_tests biomesh2_lib GTest::gtest GTest::gtest_main)
    
    # Add test
    add_test(NAME AllTests COMMAND biomesh2_tests)
else()
    message(WARNING "GoogleTest not found. Tests will not be built.")
endif()
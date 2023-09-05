cmake_minimum_required(VERSION 3.10)

set(FASTNOISE2_NOISETOOL OFF CACHE BOOL "Build Noise Tool" FORCE) #example if don't need the graph tool
add_subdirectory(deps/FastNoise2)
# find_package(FastNoise CONFIG REQUIRED) #use REQUIRED if the lib has to be found for the project to work

project(perlin_test)
add_executable(perlin_test perlin_test.cxx)

target_link_libraries(your_project PUBLIC FastNoise ...)
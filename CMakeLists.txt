cmake_minimum_required (VERSION 3.4)

project (seqan3_DNA CXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -pthread -fopenmp ")

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/submodules/seqan3/build_system")

add_library(DNA_COMMON
	common/common.h
	common/common.cpp
)
target_link_libraries (DNA_COMMON seqan3::seqan3)
target_include_directories(DNA_COMMON PRIVATE common)



add_executable (fm_indexer fm_indexer/main.cpp)
target_link_libraries (fm_indexer seqan3::seqan3)
target_include_directories(fm_indexer PRIVATE common)

add_executable (count count/main.cpp)
target_link_libraries (count seqan3::seqan3)
target_include_directories(count PRIVATE common)

add_executable (onlyhuman task/main.cpp)
target_link_libraries (onlyhuman seqan3::seqan3)

add_executable (common-breadth-serial task/Common-Breadth-Serial.cpp)
target_link_libraries (common-breadth-serial seqan3::seqan3)
target_include_directories(common-breadth-serial PRIVATE common)

add_executable (common-breadth-parallel task/Common-Breadth-Parallel.cpp)
target_link_libraries (common-breadth-parallel seqan3::seqan3)
target_include_directories(common-breadth-parallel PRIVATE common)

add_executable (temp task/temporaries.cpp)
target_link_libraries (temp seqan3::seqan3)
target_include_directories(temp PRIVATE common)
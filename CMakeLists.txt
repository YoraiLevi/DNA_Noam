cmake_minimum_required (VERSION 3.4)

project (seqan3_DNA CXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/submodules/seqan3/build_system")
 
add_executable (fm_indexer fm_indexer/main.cpp)
target_link_libraries (fm_indexer seqan3::seqan3)

add_executable (onlyhuman task/main.cpp)
target_link_libraries (onlyhuman seqan3::seqan3)

add_executable (temp task/temporaries.cpp)
target_link_libraries (temp seqan3::seqan3)
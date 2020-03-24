cmake_minimum_required (VERSION 3.4)

project (seqan3_DNA CXX)
 
find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/submodules/seqan3/build_system")
 
add_executable (fm_indexer fm_indexer/main.cpp)
target_link_libraries (fm_indexer seqan3::seqan3)

add_executable (onlyhuman task/main.cpp)
target_link_libraries (onlyhuman seqan3::seqan3)
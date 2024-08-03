cmake_minimum_required(VERSION 3.14)
project(VDP)
set(CMAKE_CXX_STANDARD 17)

find_package(TBB REQUIRED)

include_directories(src)
include_directories(vdp)
include_directories(mcl/include)

add_subdirectory(src)
add_subdirectory(vdp)
add_subdirectory(mcl)
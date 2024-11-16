aux_source_directory(. LIB_SRC)
list(FILTER LIB_SRC EXCLUDE REGEX "main-*.cpp")

add_library(vrr_lib vrr.hpp vrr.cpp)

add_library(vddlm_lib vddlm.hpp vddlm.cpp)

add_executable(vrr main-vrr.cpp)
target_link_libraries(vrr vrr_lib vdp_lib mcl mclbn384_256 TBB::tbb)

add_executable(vddlm main-vddlm.cpp)
target_link_libraries(vddlm vddlm_lib vdp_lib mcl mclbn384_256 TBB::tbb)
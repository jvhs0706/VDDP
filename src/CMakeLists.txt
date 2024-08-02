aux_source_directory(. LIB_SRC)
list(FILTER LIB_SRC EXCLUDE REGEX "main-*.cpp")

add_library(vrr_lib vrr.hpp vrr.cpp)

add_library(vdlm_lib vdlm.hpp vdlm.cpp)

add_executable(vrr main-vrr.cpp)
target_link_libraries(vrr vrr_lib vdp_lib mcl mclbn384_256)

add_executable(vdlm main-vdlm.cpp)
target_link_libraries(vdlm vdlm_lib vdp_lib mcl mclbn384_256)
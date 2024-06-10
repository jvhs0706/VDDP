aux_source_directory(. LIB_SRC)
list(FILTER LIB_SRC EXCLUDE REGEX "main*")

add_library(vdp_lib ${LIB_SRC})

add_executable(vdp_main main.cpp)
target_link_libraries(vdp_main vdp_lib mcl mclbn384_256)
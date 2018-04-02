
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/flint-fft)
include_directories(${PROJECT_SOURCE_DIR}/utils)
cado_define_test(test_transform_length.cpp LIBRARIES flint-fft utils)

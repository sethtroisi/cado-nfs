

add_library(cantor128 cantor128.c)
add_executable(cantor_bench EXCLUDE_FROM_ALL cantor_bench.c)
add_executable(cantor EXCLUDE_FROM_ALL cantor_main.c)
target_link_libraries(cantor_bench cantor128 gmp)
target_link_libraries(cantor cantor128 gmp)

add_custom_target(test-cantor
    COMMAND CANTOR=${CMAKE_CURRENT_BINARY_DIR}/cantor ${CMAKE_CURRENT_SOURCE_DIR}/check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS cantor
)

#add_library(cantor_gm128 cantor128.c mpfq_2_128.c)
#set_property(TARGET cantor_gm128 APPEND PROPERTY COMPILE_FLAGS -DCANTOR_GM)
#add_executable(cantor_gm_bench cantor_bench.c)
#add_executable(cantor_gm cantor_main.c)
#target_link_libraries(cantor_gm_bench cantor_gm128 gmp)
#target_link_libraries(cantor_gm cantor_gm128 gmp)
#
#add_custom_command(TARGET test-cantor_gm
#    COMMAND CANTOR=${CMAKE_CURRENT_BINARY_DIR}/cantor_gm ${CMAKE_CURRENT_SOURCE_DIR}/check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    DEPENDS cantor_gm
#)
#
#add_library(cantor_gmt128 cantor128.c mpfq_2_128.c)
#set_property(TARGET cantor_gmt128 APPEND PROPERTY COMPILE_FLAGS "-DCANTOR_GM -DCANTOR_GM_TRUNCATE")
#add_executable(cantor_gmt_bench cantor_bench.c)
#add_executable(cantor_gmt cantor_main.c)
#target_link_libraries(cantor_gmt_bench cantor_gmt128 gmp)
#target_link_libraries(cantor_gmt cantor_gmt128 gmp)
#
#add_custom_command(TARGET test-cantor_gmt
#    COMMAND CANTOR=${CMAKE_CURRENT_BINARY_DIR}/cantor_gmt ${CMAKE_CURRENT_SOURCE_DIR}/check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    DEPENDS cantor_gmt
#)

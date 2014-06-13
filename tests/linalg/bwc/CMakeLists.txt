
# Build makefb and las before running check
add_dependencies(check bwc_full)

macro(string_join var delimiter str_list)
    set(_ret "")
    foreach(_str ${str_list})
        if(_ret STREQUAL "")
            set(_ret "${_str}")
        else(_ret STREQUAL "")
            set(_ret "${_ret}${delimiter}${_str}")
        endif(_ret STREQUAL "")
    endforeach(_str ${str_list})
    foreach(_str ${ARGN})
        if(_ret STREQUAL "")
            set(_ret "${_str}")
        else(_ret STREQUAL "")
            set(_ret "${_ret}${delimiter}${_str}")
        endif(_ret STREQUAL "")
    endforeach(_str ${ARGN})
    set(${var} "${_ret}")
endmacro(string_join var delimiter str_list)

macro(add_bwc_test nm script)
    string_join(tail " " ${ARGN})
    add_test(test-bwc-${nm} bash -c
        "${CMAKE_CURRENT_SOURCE_DIR}/bwc-tests.sh ${CADO_NFS_BINARY_DIR}/linalg/bwc ${CADO_NFS_SOURCE_DIR}/linalg/bwc/${script} ${tail}")
endmacro(add_bwc_test)

# add_bwc_test(modp-with-sm bwc-ptrace.sh prime=65537 random_matrix_size=500 nsm=4 m=8 n=8)
# add_dependencies(check plingen_p_1)

add_bwc_test(mod2 bwc-trace.sh random_matrix_size=800 m=64 n=64)
add_dependencies(check mf_scan plingen_p_1)

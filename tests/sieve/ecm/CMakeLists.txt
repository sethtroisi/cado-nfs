# Build testbench before running check

# Get path to the testbench binary and store in variable testbench_path
# see comment in tests/sieve/CMakeLists.txt
# get_property(testbench_path TARGET testbench PROPERTY LOCATION)
set(testbench_path ${CADO_NFS_BINARY_DIR}/sieve/ecm/testbench)

# This is just to factor out the common dependency on testbench
cado_define_test(test_factor
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh

    NO_DEFAULT_RUN
    TARGET_DEPENDENCIES testbench)

cado_divert_test(test_factor pm1      ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_pm1_2_100_1000.inp -q -vf -vcf -pm1 100 1000)
cado_divert_test(test_factor pp1_27   ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_pp1_27_100_1000.inp -q -vf -vcf -pp1_27 100 1000)
cado_divert_test(test_factor pp1_65   ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_pp1_65_100_1000.inp -q -vf -vcf -pp1_65 100 1000)
cado_divert_test(test_factor ecm_10   ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecm_10_100_1000.inp -q -vf -vcf -ecm 100 1000 10)
cado_divert_test(test_factor ecm_11   ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecm_11_100_1000.inp -q -vf -vcf -ecm 100 1000 11)
cado_divert_test(test_factor ecmm12_2 ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmm12_2_100_1000.inp -q -vf -vcf -ecmm12 100 1000 2)
cado_divert_test(test_factor ecmm12_3 ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmm12_4_100_1000.inp -q -vf -vcf -ecmm12 100 1000 4)
cado_divert_test(test_factor ecmm16_1 ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmm16_1_100_1000.inp -q -vf -vcf -ecmm16 100 1000 1)
# ecmem12 with param=1 is the same as ecm with param=11
cado_divert_test(test_factor ecmem12_1 ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecm_11_100_1000.inp -q -vf -vcf -ecmem12 100 1000 1)
cado_divert_test(test_factor ecmem12_2 ${testbench_path} ${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmem12_2_100_1000.inp -q -vf -vcf -ecmem12 100 1000 2)

macro(add_order_test METHOD SIGMA PMIN PMAX)
    set(x order_ecm${METHOD}_${SIGMA}_${PMIN}_${PMAX})
    cado_divert_test(test_factor ${x}
        ${testbench_path}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_${x}.inp
        -q -vf -po${METHOD} ${SIGMA})
endmacro()

add_order_test("" 10 1000 2000)
add_order_test("" 10 1000000 1001000)
add_order_test("" 10 1000000000 1000001000)
add_order_test("" 11 1000 2000)
add_order_test("" 11 1000000 1001000)
add_order_test("" 11 1000000000 1000001000)
add_order_test("" 12 1000 2000)
add_order_test("" 12 1000000 1001000)
add_order_test("" 12 1000000000 1000001000)
add_order_test(m12 4 1000 2000)
add_order_test(m12 4 1000000 1001000)
add_order_test(m12 4 1000000000 1000001000)
add_order_test(m12 5 1000 2000)
add_order_test(m12 5 1000000 1001000)
add_order_test(m12 5 1000000000 1000001000)
add_order_test(m12 6 1000 2000)
add_order_test(m12 6 1000000 1001000)
add_order_test(m12 6 1000000000 1000001000)
add_order_test(m16 1 1000 2000)
add_order_test(m16 1 1000000 1001000)
add_order_test(m16 1 1000000000 1000001000)

add_order_test(em12 2 1000 2000)
add_order_test(em12 2 1000000 1001000)
add_order_test(em12 2 1000000000 1000001000)
add_order_test(em12 3 1000 2000)
add_order_test(em12 3 1000000 1001000)
add_order_test(em12 3 1000000000 1000001000)
add_order_test(em12 4 1000 2000)
add_order_test(em12 4 1000000 1001000)
add_order_test(em12 4 1000000000 1000001000)

cado_define_test(test_mpqs.c LIBRARIES tests facul utils)

cado_define_test(test_chains.c LIBRARIES tests facul utils)

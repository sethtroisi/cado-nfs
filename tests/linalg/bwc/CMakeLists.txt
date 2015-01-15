
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
    # Don't split this line !
    add_test(test-bwc-${nm} bash -c
        "${CMAKE_CURRENT_SOURCE_DIR}/bwc-tests.sh ${CADO_NFS_SOURCE_DIR}/linalg/bwc/${script} bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc pre_wipe=1 ${tail} -- verbose_flags=^all-bwc-sub-timings,^all-cmdline,^all-bwc-dispatch,^bwc-loading-mksol-files,^bwc-cache-major-info,^bwc-cpubinding,^perl-checks,^perl-sections")
endmacro(add_bwc_test)


# All this is 64-bit only !

set(prime_1 281474976710677)
set(prime_2 5192296858534827628530496329220121)
set(prime_3 95780971304118053647396689196894323976171195136475563)
set(prime_4 1766847064778384329583297500742918515827483896875618958121606201292619891)
set(prime_5 32592575621351777380295131014550050576823494298654980010178247189670100796213387298934358053)
set(prime_6 601226901190101306339707032778070279008174732520529886901066488712245510429339761526706943586500787976175353983)
set(prime_7 11090678776483259438313656736572334813745748301503266300681918322458485231222502492159897624416558312389564843845614287315896632389)
set(prime_8 204586912993508866875824356051724947013540127877691549342705710506008362275292159680204380770369009821930417757972504438076078534117837065833032974919)
set(prime_9 3773962424821541352241554580988268890916921220416440428376206300245624162392148852086126725177658767541468375030763844899770584629924792632561434251432696043649395327187)
set(prime_z 69617318994479297159441705409245167921344429126717528237597542082203295398081625160307507496908132931192662194421301381083506846944815643283884602656894137393981852330936660004926669193667)


macro(add_plingen_test testname)
    string_join(tail " " ${ARGN})
    # Don't split this line !
    add_test(test-bwc-plingen-${gfp_layer}-${testname} bash -c
        "${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh ${CADO_NFS_BINARY_DIR} ${tail}")
endmacro(add_plingen_test testname)

macro(add_lingen_binary_test testname)
    string_join(tail " " ${ARGN})
    # Don't split this line !
    add_test(test-bwc-lingen-binary-${testname} bash -c "${CMAKE_CURRENT_SOURCE_DIR}/test-lingen.sh ${CADO_NFS_BINARY_DIR} ${tail}")
endmacro(add_lingen_binary_test)

add_lingen_binary_test(test0 1dea0d8879d125aa280e26b8573f932f6bccea75 64 64 160 1 lingen-threshold=10 ${extra})

foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(extra)
    string(REGEX MATCH "p_([0-9]+)" t "${gfp_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    else(t)
        if(gfp_layer STREQUAL "pz")
            set(t TRUE)
            set(width "z")
            set(extra "plingen_program=plingen_pz")
        endif(gfp_layer STREQUAL "pz")
    endif(t)
    set(p "${prime_${width}}")
    add_bwc_test(modp-inhomogeneous-${gfp_layer} bwc-ptrace.sh prime=${p} random_matrix_size=500 nrhs=4 m=8 n=8 seed=1 ${extra})
    if(gfp_layer STREQUAL "p_1")
        set(mextra)
        # Only try to enable magma if CHECKS_EXPENSIVE is set.
        if(DEFINED ENV{CHECKS_EXPENSIVE} AND DEFINED ENV{MAGMA})
            set(mextra magma=$ENV{MAGMA})
        endif(DEFINED ENV{CHECKS_EXPENSIVE} AND DEFINED ENV{MAGMA})
        add_bwc_test(modp-inhomogeneous-minimal-3rhs-mn6 bwc-ptrace.sh prime=1009
            random_matrix_size=50 nrhs=3 m=6 n=6 seed=1 ${mextra})
        add_bwc_test(modp-inhomogeneous-minimal-2rhs-mn4 bwc-ptrace.sh
            prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1 ${mextra})
        if(DEFINED ENV{CHECKS_EXPENSIVE})
            add_bwc_test(modp-inhomogeneous-minimal-2rhs-mn4-expensive bwc-ptrace.sh
                prime=1237940039285380274899124357 random_matrix_size=200 nrhs=2 m=4 n=4 seed=1 ${mextra})
        endif(DEFINED ENV{CHECKS_EXPENSIVE})

        add_bwc_test(modp-inhomogeneous-minimal-1rhs-mn1 bwc-ptrace.sh
            prime=1009 random_matrix_size=50 nrhs=1 m=1 n=1 seed=1 ${mextra})
#        add_bwc_test(modp-homogeneous-minimal-mn4 bwc-ptrace.sh
#            prime=1009 random_matrix_size=50 m=4 n=4 seed=1 ${mextra})
#        add_bwc_test(modp-inhomogeneous-minimal-2rhs-mn2 bwc-ptrace.sh
#            prime=1009 random_matrix_size=50 nrhs=2 m=2 n=2 seed=1 ${mextra})
    endif(gfp_layer STREQUAL "p_1")
    add_dependencies(check plingen_${gfp_layer})

    # Also add some plingen tests if it so occurs that the current
    # backend can process them
    # Each of the add_plingen_test lines has the following arguments
    # plingen test name
    # reference sha1 sum for generator
    # m
    # n
    # length of sequence to compute
    # p
    # seed for random generation (entries are (2^k*(seed/1000) mod q) + (3^k*(seed%1000) mod q) mod p, with q=2^20-3).
    # arguments to be passed to plingen (including mpi arguments)

    if(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_1")
        add_plingen_test(test0
            2b6428178054fc0d3931dcdc9257ccf3bdea62f9 4 2 36 1009 17042
            lingen-threshold=10 lingen-mpi-threshold=10 thr=2x2 --ascii
            ${extra})
        add_plingen_test(test1
            f9bc59f0caaca0fd523ef6a020a6f49fa97f38e9 4 2 200 1009 88888
            lingen-threshold=10 lingen-mpi-threshold=10 thr=2x2 --ascii
            ${extra})
        add_plingen_test(test2
            f9bc59f0caaca0fd523ef6a020a6f49fa97f38e9 4 2 200 1009 88888
            lingen-threshold=10 lingen-mpi-threshold=40 thr=2x2 --ascii
            ${extra})
        add_plingen_test(test3
            b448493e04299d8d54be6cc1273c43293520ac9e 3 3 1000 1009 1111
            lingen-threshold=10 lingen-mpi-threshold=40 thr=3x3 --ascii
            ${extra})
    endif(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_1")
    if(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_3")
        add_plingen_test(test4
            00a481d0e0e2ee7403fa234b7f8b818dc2ec8db4 3 3 1000
            1532495540865888858358347027150309183618739122183602191 1111
            lingen-threshold=10 lingen-mpi-threshold=40 thr=3x3
            ${extra})
    endif(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_3")
endforeach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})


add_bwc_test(mod2 bwc-trace.sh random_matrix_size=800 m=64 n=64 seed=1 thr=2x2)
add_dependencies(check mf_scan)

add_bwc_test(mod2-oddthr bwc-trace.sh random_matrix_size=800 mn=128 seed=1 thr=3)

add_bwc_test(mod2-interleaving bwc-trace.sh random_matrix_size=800 mn=128
    interleaving=1 seed=1 thr=2)

add_bwc_test(mod2-m-neq-n bwc-trace.sh random_matrix_size=400 m=128 n=64
    seed=1 thr=5)

add_subdirectory(mpfq)

if(HAVE_CXX11 AND HAVE_HWLOC)
    include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc)
    build_and_run_test(test_cpubinding.c "utils;bwc_base" --cpubinding ${CMAKE_CURRENT_SOURCE_DIR}/cpubind_config_file_for_tests.conf)
endif(HAVE_CXX11 AND HAVE_HWLOC)

macro(test_dispatch_matmul testname)
    string_join(tail " " ${ARGN})
    # Don't split this line !
    add_test(${testname} bash -c "${CMAKE_CURRENT_SOURCE_DIR}/bwc-dispatch-matmul.sh bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc ${tail}")
endmacro(test_dispatch_matmul testname)

add_dependencies(check bwc_matmul)
# mingw has antique bash (3.1) with different regexp syntax, and having
# the script below work with both isn't a very funny perspective.
if (NOT MINGW)
test_dispatch_matmul(dispatch-matmul-65536-thr nrows=65536 density=2 thr=2x2)
if (HAVE_MPI)
test_dispatch_matmul(dispatch-matmul-65536-mpi nrows=65536 density=2 mpi=2x2)
endif (HAVE_MPI)
endif (NOT MINGW)


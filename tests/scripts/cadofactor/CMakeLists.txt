

set(PYTHON_DOCTEST_SOURCES
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadoparams.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadoprograms.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadotask.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/workunit.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/wudb.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/wuserver.py
)

# make sure we byte compile everything before running the doctests.
cado_define_test(TEST_NAME python_compile
    PROGRAM python3 -m py_compile ${PYTHON_DOCTEST_SOURCES})

foreach (doctest_source ${PYTHON_DOCTEST_SOURCES})
  # Get filename without extension, e.g., cadoparams
  get_filename_component (basename ${doctest_source} NAME_WE)
  # Add doctest of Python sources, e.g., doctest cadoparams.py under test
  # name cadoparams
  cado_define_test(TEST_NAME test_python_${basename}
      PROGRAM python3 -m doctest ${doctest_source}
      TEST_DEPENDENCIES python_compile)
endforeach (doctest_source)

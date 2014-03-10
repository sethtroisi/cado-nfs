FILE(GLOB pythonscripts "${CMAKE_CURRENT_SOURCE_DIR}/cado*.py" "${CMAKE_CURRENT_SOURCE_DIR}/wu*.py")
install(PROGRAMS ${pythonscripts} patterns.py workunit.py upload.py DESTINATION bin)

set(PYTHON_DOCTEST_SOURCES
  cadoparams.py
  cadoprograms.py
  cadotask.py
  workunit.py
  wudb.py
  wuserver.py
)

foreach (doctest_source ${PYTHON_DOCTEST_SOURCES})
  # Get filename without extension, e.g., cadoparams
  get_filename_component (basename ${doctest_source} NAME_WE)
  # Add doctest of Python sources, e.g., doctest cadoparams.py under test
  # name cadoparams
  add_test(${basename} python3 -m doctest ${CMAKE_CURRENT_SOURCE_DIR}/${doctest_source})
endforeach (doctest_source)

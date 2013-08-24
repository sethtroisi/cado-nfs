FILE(GLOB pythonscripts "${CMAKE_CURRENT_SOURCE_DIR}/cado*.py" "${CMAKE_CURRENT_SOURCE_DIR}/wu*.py")
install(PROGRAMS ${pythonscripts} patterns.py workunit.py upload.py DESTINATION bin)

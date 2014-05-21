
configure_file(cadofactor.py cadofactor.py ESCAPE_QUOTES @ONLY)
configure_file(wudb.py wudb.py ESCAPE_QUOTES @ONLY)
configure_file(wuclient2.py wuclient2.py ESCAPE_QUOTES @ONLY)

install(FILES cadocommand.py DESTINATION ${LIBSUFFIX})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cadofactor.py DESTINATION ${BINSUFFIX})
install(FILES cadologger.py DESTINATION ${LIBSUFFIX})
install(FILES cadoparams.py DESTINATION ${LIBSUFFIX})
install(FILES cadoprograms.py DESTINATION ${LIBSUFFIX})
install(FILES cadotask.py DESTINATION ${LIBSUFFIX})
install(FILES patterns.py DESTINATION ${LIBSUFFIX})
# install(FILES signalhandler.py DESTINATION ${LIBSUFFIX})
install(FILES upload.py DESTINATION ${LIBSUFFIX})
# install(FILES verify.py DESTINATION ${LIBSUFFIX})
install(FILES workunit.py DESTINATION ${LIBSUFFIX})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/wuclient2.py DESTINATION ${BINSUFFIX})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/wudb.py DESTINATION ${BINSUFFIX})
install(FILES wuserver.py DESTINATION ${LIBSUFFIX})

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

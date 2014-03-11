# Build sm before running check
add_dependencies(check sm)

# Get path to the sm binary and store in variable sm_path
get_property(sm_path TARGET sm PROPERTY LOCATION)

add_test(testsm bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/testsm.sh\" \"${sm_path}\"")

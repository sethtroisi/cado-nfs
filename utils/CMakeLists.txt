
SET (sources 
    polyfile.c
    relation.c
    timing.c
    long_poly.c
    getprime.c
    mod_ul.c)

ADD_LIBRARY (utils ${sources})

ADD_EXECUTABLE (cache cache.c)
TARGET_LINK_LIBRARIES(cache utils)


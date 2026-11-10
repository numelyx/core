if(NOT DEFINED PREFIX)
    message(FATAL_ERROR "PREFIX must point to the installed package root")
endif()

foreach(header IN ITEMS
    numelyx.h
    export.h
    version.h
    status.h
    matrix.h
    rng.h
    optimizer.h
    probability.h
)
    if(NOT EXISTS "${PREFIX}/include/numelyx/${header}")
        message(FATAL_ERROR "Missing installed public header: ${header}")
    endif()
endforeach()

if(EXISTS "${PREFIX}/include/numelyx/detail")
    message(FATAL_ERROR "Private implementation headers must not be installed")
endif()

if(NOT DEFINED NM OR NM STREQUAL "")
    message(FATAL_ERROR "NM must name the platform symbol inspection tool")
endif()
if(NOT DEFINED LIBRARY OR NOT EXISTS "${LIBRARY}")
    message(FATAL_ERROR "LIBRARY must name an existing shared library")
endif()

execute_process(
    COMMAND "${NM}" -D --defined-only "${LIBRARY}"
    RESULT_VARIABLE nm_status
    OUTPUT_VARIABLE symbols
    ERROR_VARIABLE nm_error
)
if(NOT nm_status EQUAL 0)
    message(FATAL_ERROR "Could not inspect shared-library symbols: ${nm_error}")
endif()

if(NOT symbols MATCHES "[ \t]numelyx_matrix_create_dense(@@[^ \t\n]+)?([ \t\n]|$)")
    message(FATAL_ERROR "Supported public symbol numelyx_matrix_create_dense is not exported")
endif()

if(symbols MATCHES "[ \t]numelyx_core_matrix_create(@@[^ \t\n]+)?([ \t\n]|$)")
    message(FATAL_ERROR "Private Core symbol numelyx_core_matrix_create leaked from shared library")
endif()

cmake_minimum_required(VERSION 3.13)
project(ts LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(
        CMAKE_BUILD_TYPE
        "Release"
        CACHE STRING "Choose the build type."
        FORCE
    )
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(Threads REQUIRED)

# -----------------------------------------------------------------------------
# Gurobi 12.0
# -----------------------------------------------------------------------------

if(NOT DEFINED ENV{GUROBI_HOME} OR "$ENV{GUROBI_HOME}" STREQUAL "")
    message(FATAL_ERROR
        "GUROBI_HOME is not set. Point it to Gurobi's linux64 directory."
    )
endif()

file(TO_CMAKE_PATH "$ENV{GUROBI_HOME}" GUROBI_HOME)

find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c++.h
    PATHS "${GUROBI_HOME}/include"
    NO_DEFAULT_PATH
)

find_library(GUROBI_LIBRARY
    NAMES gurobi120
    PATHS "${GUROBI_HOME}/lib"
    NO_DEFAULT_PATH
)

# Prefer Gurobi's generic C++ library name. The second name supports the
# compiler-specific archive used by the original project configuration.
find_library(GUROBI_CXX_LIBRARY
    NAMES gurobi_c++ gurobi_g++8.5
    PATHS "${GUROBI_HOME}/lib"
    NO_DEFAULT_PATH
)

if(NOT GUROBI_INCLUDE_DIR)
    message(FATAL_ERROR
        "Could not find gurobi_c++.h under "
        "GUROBI_HOME='${GUROBI_HOME}'."
    )
endif()

if(NOT GUROBI_LIBRARY)
    message(FATAL_ERROR
        "Could not find the Gurobi 12.0 library under "
        "GUROBI_HOME='${GUROBI_HOME}'."
    )
endif()

if(NOT GUROBI_CXX_LIBRARY)
    message(FATAL_ERROR
        "Could not find the Gurobi C++ library under "
        "GUROBI_HOME='${GUROBI_HOME}'."
    )
endif()

message(STATUS "Gurobi home: ${GUROBI_HOME}")
message(STATUS "Gurobi include directory: ${GUROBI_INCLUDE_DIR}")
message(STATUS "Gurobi library: ${GUROBI_LIBRARY}")
message(STATUS "Gurobi C++ library: ${GUROBI_CXX_LIBRARY}")

# -----------------------------------------------------------------------------
# LEMON 1.3.1
# -----------------------------------------------------------------------------

# Supply these paths when configuring:
#
#   -DLEMON_SOURCE_DIR=/path/to/lemon-1.3.1
#   -DLEMON_BUILD_DIR=/path/to/lemon-1.3.1/build

set(LEMON_SOURCE_DIR
    ""
    CACHE PATH "LEMON source directory"
)

set(LEMON_BUILD_DIR
    ""
    CACHE PATH "LEMON build directory"
)

if(LEMON_SOURCE_DIR STREQUAL "")
    message(FATAL_ERROR
        "LEMON_SOURCE_DIR is not set. Configure the project with "
        "-DLEMON_SOURCE_DIR=/path/to/lemon-1.3.1."
    )
endif()

# If only the source directory is provided, use its build/ subdirectory.
if(LEMON_BUILD_DIR STREQUAL "")
    set(
        LEMON_BUILD_DIR
        "${LEMON_SOURCE_DIR}/build"
        CACHE PATH "LEMON build directory"
        FORCE
    )
endif()

# The source tree contains ordinary headers.
if(NOT EXISTS "${LEMON_SOURCE_DIR}/lemon/list_graph.h")
    message(FATAL_ERROR
        "LEMON source headers were not found under "
        "'${LEMON_SOURCE_DIR}'."
    )
endif()

# The build tree contains generated headers such as lemon/config.h.
if(NOT EXISTS "${LEMON_BUILD_DIR}/lemon/config.h")
    message(FATAL_ERROR
        "LEMON generated headers were not found under "
        "'${LEMON_BUILD_DIR}'. Build LEMON first or correct "
        "LEMON_BUILD_DIR."
    )
endif()

# Upstream LEMON 1.3.1 normally produces lemon/libemon.a. The alternative
# name supports distributions that rename the library to liblemon.
find_library(LEMON_LIBRARY
    NAMES emon lemon
    PATHS
        "${LEMON_BUILD_DIR}/lemon"
        "${LEMON_BUILD_DIR}/lib"
        "${LEMON_BUILD_DIR}"
    NO_DEFAULT_PATH
)

if(NOT LEMON_LIBRARY)
    message(FATAL_ERROR
        "Could not find the compiled LEMON library under "
        "'${LEMON_BUILD_DIR}'. Build LEMON first or correct "
        "LEMON_BUILD_DIR."
    )
endif()

get_filename_component(
    LEMON_LIBRARY_DIR
    "${LEMON_LIBRARY}"
    DIRECTORY
)

message(STATUS "LEMON source directory: ${LEMON_SOURCE_DIR}")
message(STATUS "LEMON build directory: ${LEMON_BUILD_DIR}")
message(STATUS "LEMON library: ${LEMON_LIBRARY}")

# -----------------------------------------------------------------------------
# Project sources
# -----------------------------------------------------------------------------

# src/main.cpp contains main_mink() and main_maxp(). The two wrapper files
# provide the standard main() entry point for their respective executables.
set(
    TS_DRIVER_SOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)

set(
    TS_MINK_ENTRY_SOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main_mink_entry.cpp"
)

set(
    TS_MAXP_ENTRY_SOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main_maxp_entry.cpp"
)

foreach(REQUIRED_SOURCE
        "${TS_DRIVER_SOURCE}"
        "${TS_MINK_ENTRY_SOURCE}"
        "${TS_MAXP_ENTRY_SOURCE}")
    if(NOT EXISTS "${REQUIRED_SOURCE}")
        message(FATAL_ERROR
            "Required source file not found: ${REQUIRED_SOURCE}"
        )
    endif()
endforeach()

# CONFIGURE_DEPENDS causes CMake to detect added or removed .cpp files when
# the generated build system is invoked.
file(GLOB_RECURSE TS_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# These files are added explicitly to the appropriate executable targets.
# Removing them from the shared list prevents duplicate definitions.
list(REMOVE_ITEM TS_SOURCES
    "${TS_DRIVER_SOURCE}"
    "${TS_MINK_ENTRY_SOURCE}"
    "${TS_MAXP_ENTRY_SOURCE}"
)

add_executable(ts_mink
    ${TS_SOURCES}
    "${TS_DRIVER_SOURCE}"
    "${TS_MINK_ENTRY_SOURCE}"
)

add_executable(ts_maxp
    ${TS_SOURCES}
    "${TS_DRIVER_SOURCE}"
    "${TS_MAXP_ENTRY_SOURCE}"
)

# -----------------------------------------------------------------------------
# Shared target configuration
# -----------------------------------------------------------------------------

function(configure_ts_target TARGET_NAME)
    target_include_directories(${TARGET_NAME} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/common"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/single_telescope"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/multi_telescope/min_telescopes"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/multi_telescope/max_probabilities"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/apps"
        "${GUROBI_INCLUDE_DIR}"
        "${LEMON_SOURCE_DIR}"
        "${LEMON_BUILD_DIR}"
    )

    # The Gurobi C++ wrapper is normally a static archive that refers to
    # symbols in the main Gurobi library, so it is listed first.
    target_link_libraries(${TARGET_NAME} PRIVATE
        OpenMP::OpenMP_CXX
        "${GUROBI_CXX_LIBRARY}"
        "${GUROBI_LIBRARY}"
        "${LEMON_LIBRARY}"
        Threads::Threads
        "${CMAKE_DL_LIBS}"
        m
    )

    # GCC 8 and older require a separate library for std::filesystem.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
       AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
        target_link_libraries(${TARGET_NAME} PRIVATE stdc++fs)
    endif()

    # Allow the binaries to locate shared dependency libraries without
    # modifying LD_LIBRARY_PATH.
    set_target_properties(${TARGET_NAME} PROPERTIES
        BUILD_RPATH
            "${GUROBI_HOME}/lib;${LEMON_LIBRARY_DIR}"
        INSTALL_RPATH
            "${GUROBI_HOME}/lib;${LEMON_LIBRARY_DIR}"
    )
endfunction()

configure_ts_target(ts_mink)
configure_ts_target(ts_maxp)
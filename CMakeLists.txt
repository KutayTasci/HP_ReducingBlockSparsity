cmake_minimum_required(VERSION 3.10)
project(HypergraphPartitioner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------
# AddressSanitizer option
# -------------------------------
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

if(ENABLE_ASAN)
    message(STATUS "Building with AddressSanitizer")
    set(SANITIZER_FLAGS "-fsanitize=address -fno-omit-frame-pointer -O1 -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

# -------------------------------
# Detect OpenMP
# -------------------------------
find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()

# -------------------------------
# Source files
# -------------------------------
file(GLOB HPNBM_SOURCES "${CMAKE_SOURCE_DIR}/src/HPNBM/*.cpp")

set(SOURCES
    main.cc
    ${HPNBM_SOURCES}
)

# -------------------------------
# Include directories
# -------------------------------
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    /opt/boost_1_82_0/include
    ./ext/      # add PaToH headers
)

link_directories(
    /opt/boost_1_82_0/lib
    /usr/local/lib
    ./ext/
)

# -------------------------------
# Executable
# -------------------------------
add_executable(TestHP ${SOURCES})

target_link_libraries(TestHP 
    mtkahypar
    boost_program_options
    OpenMP::OpenMP_CXX   # OpenMP
    patoh                # PaToH
)

# -------------------------------
# RPATH settings
# -------------------------------
set_target_properties(TestHP PROPERTIES
    BUILD_RPATH "/opt/boost_1_82_0/lib;/usr/local/lib;./ext/"
    INSTALL_RPATH "/opt/boost_1_82_0/lib;/usr/local/lib;./ext/"
)

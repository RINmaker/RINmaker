cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

# require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "you cannot build in a source directory (or any directory with a CMakeLists.txt file). please make a build subdirectory. feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

project(RINmaker)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# pugixml is available under MIT license (see repo)
message("fetching pugixml...")
FetchContent_Declare(
        pugixml
        GIT_REPOSITORY https://github.com/RINmaker/pugixml.git
        GIT_TAG v1.12.1
)
FetchContent_MakeAvailable(pugixml)

# CLI11 is available under their own license (see repo)
# the license permits usage and redistribution
message("fetching CLI11...")
FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/RINmaker/CLI11.git
        GIT_TAG v2.1.2
)
FetchContent_MakeAvailable(CLI11)

# spdlog is available under MIT license (see repo)
message("fetching spdlog...")
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/RINmaker/spdlog.git
        GIT_TAG v1.9.2
)
FetchContent_MakeAvailable(spdlog)

option(USE_DEBUG "enter debug mode for RINmaker_library" OFF)
if (USE_DEBUG)
    add_definitions(-DDEBUG)
endif()

add_subdirectory(app)
add_subdirectory(test)

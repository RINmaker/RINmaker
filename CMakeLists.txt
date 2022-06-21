cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 17)

project("RINmaker")

file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
            "Out-of-source builds are REQUIRED. Feel free to remove any generated `CMakeCache.txt` and `CMakeFiles`")
endif ()

include(FetchContent)
function("RINmaker_FetchContent" NAME TAG)
    message("Fetching ${NAME}...")
    FetchContent_Declare(${NAME} GIT_REPOSITORY "https://github.com/RINmaker/${NAME}.git" GIT_TAG ${TAG})
    FetchContent_MakeAvailable(${NAME})
    message("Fetch ${NAME} OK.")
endfunction()

RINmaker_FetchContent("pugixml" "v1.12.1")
RINmaker_FetchContent("CLI11" "v2.1.2")
RINmaker_FetchContent("spdlog" "v1.9.2")

add_subdirectory(app)
add_subdirectory(test)

cmake_minimum_required(VERSION 3.21)

project(RINmaker)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# pugixml is available under MIT license (see repo)
message("fetching pugixml...")
FetchContent_Declare(
        pugixml
        GIT_REPOSITORY https://github.com/zeux/pugixml.git
        GIT_TAG v1.12.1
)
FetchContent_MakeAvailable(pugixml)

# CLI11 is available under their own license (see repo)
# the license permits usage and redistribution
message("fetching CLI11...")
FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v2.1.2
)
FetchContent_MakeAvailable(CLI11)

# spdlog is available under MIT license (see repo)
message("fetching spdlog...")
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.9.2
)
FetchContent_MakeAvailable(spdlog)

# googletest is available under their own license (see repo)
# the license permits usage and redistribution
message("fetching googletest...")
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.11.0
)
FetchContent_MakeAvailable(googletest)

add_subdirectory(src)
# googletest is available under their own license (see repo)
# the license permits usage and redistribution
message("fetching googletest...")
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/RINmaker/googletest.git
        GIT_TAG release-1.11.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

file(GLOB SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB APP_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../app/sources/*.cpp)
add_executable(RINmaker_test ${SRCS} ${APP_SRCS})

target_include_directories(RINmaker_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../app/include")
target_link_libraries(RINmaker_test PRIVATE gtest spdlog pugixml CLI11)

# copy test_case folder near target executable
add_custom_command(
        TARGET RINmaker_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../test_case/ $<TARGET_FILE_DIR:RINmaker_test>/test_case/)
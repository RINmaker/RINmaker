set(TARGET_NAME "RINmaker_test")

message("Fetching googletest...")
FetchContent_Declare(
        "googletest"
        GIT_REPOSITORY "https://github.com/RINmaker/googletest.git"
        GIT_TAG "release-1.11.0")

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
message("Fetch googletest OK.")

file(GLOB SRCS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB APP_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/../app/sources/*.cpp")

add_executable(${TARGET_NAME} ${SRCS} ${APP_SRCS})

target_include_directories(
        ${TARGET_NAME}
        PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../app/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../gemmi/include")

target_link_libraries(${TARGET_NAME} PRIVATE "gtest" "spdlog" "pugixml" "CLI11")

add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/../test_case/"
        "$<TARGET_FILE_DIR:RINmaker_test>/test_case/")
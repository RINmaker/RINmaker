set(TARGET_NAME "RINmaker")

file(GLOB SRCS "${CMAKE_CURRENT_SOURCE_DIR}/sources/*.cpp")
add_executable(${TARGET_NAME} "main.cpp" ${SRCS})

target_link_libraries(${TARGET_NAME} PRIVATE "spdlog" "pugixml" "CLI11" -static)

target_include_directories(
        ${TARGET_NAME} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../gemmi/include")
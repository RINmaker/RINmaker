file(GLOB SRCS ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.cpp)
add_executable(RINmaker main.cpp ${SRCS})

target_link_libraries(RINmaker PRIVATE spdlog pugixml CLI11 -static)

target_include_directories(RINmaker PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
add_executable(RINmaker main.cpp)

target_link_libraries(RINmaker PRIVATE RINmaker_library spdlog pugixml CLI11)
target_include_directories(RINmaker PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../headers")

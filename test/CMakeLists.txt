add_executable(RINmaker_test all_tests.cpp BlackBoxTest.cpp KDTreeTest.cpp)

target_include_directories(RINmaker_test PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../headers")
target_link_libraries(RINmaker_test PRIVATE RINmaker_library gtest spdlog pugixml CLI11)

# copy test_case folder near target executable
add_custom_command(
        TARGET RINmaker_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../test_case/ $<TARGET_FILE_DIR:RINmaker_test>/test_case/)
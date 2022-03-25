set(RINMAKER_SOURCES
        computation.cpp
        config.cpp
        energy.cpp
        chemical_entities.cpp
        log_manager.cpp
        main.cpp
        bond_network.cpp
        runtime_params.cpp
        pdb_records.cpp
        pdb_data.cpp
        graphml_output.cpp
        secondary_structures.cpp)

set(TEST_SOURCES "tests/BlackBoxTest.cpp" "tests/KDTreeTest.cpp")

add_executable(RINmaker ${RINMAKER_SOURCES} ${TEST_SOURCES})
target_link_libraries(RINmaker PRIVATE spdlog pugixml gtest CLI11)
target_include_directories(RINmaker PRIVATE "../headers")

add_custom_command(
        TARGET RINmaker POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../test_case/ $<TARGET_FILE_DIR:RINmaker>/testcase/)
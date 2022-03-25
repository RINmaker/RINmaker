set(RINMAKER_HEADERS "../headers")
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

set(RINMAKER_DEPS spdlog pugixml CLI11)

set(RINMAKER_TEST_SOURCES "tests/BlackBoxTest.cpp" "tests/KDTreeTest.cpp")
set(RINMAKER_TEST_DEPS gtest)

add_executable(RINmaker ${RINMAKER_SOURCES} ${RINMAKER_TEST_SOURCES})
target_link_libraries(RINmaker PRIVATE ${RINMAKER_DEPS} ${RINMAKER_TEST_DEPS})

target_include_directories(RINmaker PRIVATE ${RINMAKER_HEADERS})

add_custom_command(
        TARGET RINmaker POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../test_case/ $<TARGET_FILE_DIR:RINmaker>/testcase/)
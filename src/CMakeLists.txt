set(APP_SRCS
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

set(TEST_SRCS "tests/BlackBoxTest.cpp" "tests/KDTreeTest.cpp")

add_executable(RINmaker ${APP_SRCS} ${TEST_SRCS})
target_link_libraries(RINmaker PRIVATE spdlog pugixml gtest CLI11)
target_include_directories(RINmaker PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/headers")

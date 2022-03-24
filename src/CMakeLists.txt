set(SOURCES
        computation.cpp
        config.cpp
        energy.cpp
        chemical_entities.cpp
        log_manager.cpp
        main.cpp
        bond_network.cpp
        runtime_params.cpp
        pdb_records.cpp
        graphml_output.cpp
        secondary_structures.cpp)

add_executable(RINmaker ${SOURCES})
target_link_libraries(RINmaker PRIVATE spdlog pugixml gtest CLI11)
target_include_directories(RINmaker PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/headers")

add_library(
        RINmaker_library

        rin_maker.cpp
        config.cpp
        energy.cpp
        chemical_entity.cpp
        log_manager.cpp
        bond_network.cpp
        runtime_params.cpp
        pdb_records.cpp
        pdb_data.cpp
        graphml_output.cpp
        secondary_structures.cpp)

target_link_libraries(RINmaker_library PRIVATE spdlog pugixml CLI11)
target_include_directories(RINmaker_library PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../headers")

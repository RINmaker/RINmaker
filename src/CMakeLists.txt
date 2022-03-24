set(SOURCES
        computation.cpp
        config.cpp
        energy.cpp
        entity.cpp
        log_manager.cpp
        main.cpp
        network.cpp
        parameters.cpp
        records.cpp
        rin.cpp
        structure.cpp)

add_executable(RINmaker ${SOURCES})
target_link_libraries(RINmaker PRIVATE spdlog pugixml gtest CLI11)
target_include_directories(RINmaker PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/headers")

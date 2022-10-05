#pragma once

#pragma warning(push, 0)

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>

#pragma warning(pop)

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>

#include <filesystem>

class log_manager
{
private:
    std::shared_ptr<spdlog::logger> main_logger{};

    log_manager() = default;

    static log_manager& instance();

public:
    static void initialize(std::filesystem::path const& log_file, bool verbose);

    static std::shared_ptr<spdlog::logger> main()
    { return instance().main_logger; }
};

#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>

#include <filesystem>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>

using std::string;

class log_manager {
private:
    //std::unordered_map<string, std::shared_ptr<spdlog::logger>> loggers;
    std::shared_ptr<spdlog::logger> console_logger;
    std::shared_ptr<spdlog::logger> file_logger;
    std::shared_ptr<spdlog::logger> main_logger;

    //std::filesystem::path log_directory;

    log_manager() = default;

    static log_manager &instance() {
        static log_manager lm;
        return lm;
    }

public:
    static void initialize(std::filesystem::path const &log_directory);   // throws

    //static void insert(string const& id, string const& sink_name);

    //static std::shared_ptr<spdlog::logger> get_default();   // TODO attualmente get_default() � cout: farei in modo che punti invece ad un logger UNICO che � quello di default
    //static std::shared_ptr<spdlog::logger> get(string const& key);
    //static std::vector<string> get_ids();

    static void mirror_in_stdout(bool b);

    static std::shared_ptr<spdlog::logger> main() { return instance().main_logger; }

    static std::shared_ptr<spdlog::logger> file() { return instance().file_logger; }

    static std::shared_ptr<spdlog::logger> console() { return instance().console_logger; }
};

// #define LOG (*log_manager::main())
